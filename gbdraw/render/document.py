"""SVG document text: the one serialization of a drawing per request."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, TextIO, TypeAlias

from svgwrite import Drawing
from svgwrite.base import BaseElement

from gbdraw.core.figure_length import FigureLength, RequestedFigureSize
from gbdraw.render.group_content import (
    GroupContent,
    element_text,
    placed_content,
    require_element_markup,
    start_tag,
)


def parse_svg_number(value: object, *, default: float = 0.0) -> float:
    if value is None:
        return float(default)
    raw = str(value).strip()
    if raw.endswith("px"):
        raw = raw[:-2]
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default)


def drawing_size(drawing: Drawing) -> tuple[float, float]:
    """A drawing's width and height in SVG px; 0 for a size that is not a number."""

    return (
        parse_svg_number(drawing.attribs.get("width")),
        parse_svg_number(drawing.attribs.get("height")),
    )


# What a drawing drawn to a requested size was asked for, so its written file states the size in that unit.
_EXPORT_SIZE_ATTR = "_gbdraw_export_size"
# The decimals a side converted to the other side's unit is written with.
_CONVERTED_DECIMALS = 4


def bind_export_size(drawing: Drawing, size: RequestedFigureSize) -> None:
    """Keep the size a drawing was fitted to, for ``write_svg_document``; the drawing itself stays in px."""

    setattr(drawing, _EXPORT_SIZE_ATTR, size if size.requested else None)


def _export_root_size(drawing: Drawing) -> dict[str, str] | None:
    """The root ``width`` and ``height`` of a written file: each side in its requested unit, else in the other's.

    The ``viewBox`` stays in px, so a viewer scales the px drawing to the physical size.
    """

    size = getattr(drawing, _EXPORT_SIZE_ATTR, None)
    if not isinstance(size, RequestedFigureSize) or not size.requested:
        return None
    unit = size.width if size.width is not None else size.height
    assert unit is not None
    width_px, height_px = drawing_size(drawing)

    def side(requested: FigureLength | None, px: float) -> str:
        if requested is not None:
            return requested.text
        converted = unit.of_px(px)
        return FigureLength(round(converted.value, _CONVERTED_DECIMALS), converted.unit).text

    return {"width": side(size.width, width_px), "height": side(size.height, height_px)}


def _root_with(drawing: Drawing, attribs: dict[str, str] | None) -> Drawing:
    """``drawing`` itself, or a shallow copy whose root carries ``attribs``; the drawing is not changed."""

    if not attribs:
        return drawing
    root = object.__new__(type(drawing))
    root.__dict__.update(drawing.__dict__)
    root.attribs = {**drawing.attribs, **attribs}
    return root


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


def document_parts(drawing: Drawing, *, root_attribs: dict[str, str] | None = None) -> DocumentParts:
    """Serialize ``drawing`` once, as ordered pieces (the one order of serialization).

    ``root_attribs`` replace root attributes in the text only (the written size of ``write_svg_document``).

    The prolog is svgwrite's: the XML declaration and the stylesheet
    instructions. A placed content group with children is a ``ContentPiece``
    over the stored text of its children; every other top-level element is a
    ``RawPiece`` by svgwrite's own serialization. Raises ``ExportError`` when a
    piece holds a comment, CDATA or a processing instruction.
    """

    prolog = _DECLARATION + "".join(
        _STYLESHEET % stylesheet for stylesheet in drawing._stylesheets
    )
    root = _root_with(drawing, root_attribs)
    if not drawing.elements:
        return DocumentParts(
            prolog=prolog, head=require_element_markup(root.tostring()), pieces=(), tail=""
        )
    return DocumentParts(
        prolog=prolog,
        head=require_element_markup(start_tag(root)),
        pieces=tuple(_piece(element) for element in drawing.elements),
        tail=f"</{drawing.elementname}>",
    )


def write_svg_document(drawing: Drawing) -> SvgDocument:
    """Serialize ``drawing`` once, as ``Drawing.write`` would write it (``document_parts``).

    A drawing fitted to a requested size (``bind_export_size``) states it in the requested unit.
    """

    parts = document_parts(drawing, root_attribs=_export_root_size(drawing))
    return SvgDocument(prolog=parts.prolog, body=parts.body())


__all__ = [
    "ContentPiece",
    "DocumentParts",
    "Piece",
    "RawPiece",
    "SvgDocument",
    "bind_export_size",
    "document_parts",
    "write_svg_document",
    "write_text",
]
