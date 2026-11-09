"""SVG document text: the one serialization of a drawing per request."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
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


# A requested side the drawn canvas kept: within this many px it is written as requested.
_KEPT_SIDE_PX = 0.01


def export_root_size(size: RequestedFigureSize, width_px: float, height_px: float) -> dict[str, str] | None:
    """The root ``width`` and ``height`` a written file of a ``width_px`` x ``height_px`` drawing fitted to ``size``
    states; ``None`` for a size with no requested side (the root stays in px).

    A requested side the canvas kept is written as requested; any other side is converted, into its own requested
    unit (a canvas that grew past the box after the fit) or into the first requested unit, and rounded to
    ``_CONVERTED_DECIMALS``. The ``viewBox`` stays in px, so a viewer scales the px drawing to the physical size.
    The one owner for the CLI and Python export and the Web SVG export (``gbdraw/web_support/export_limit.py``).
    """

    if not size.requested:
        return None
    first = size.width if size.width is not None else size.height
    assert first is not None

    def side(requested: FigureLength | None, px: float) -> str:
        if requested is not None and abs(requested.px - px) <= _KEPT_SIDE_PX:
            return requested.text
        converted = (requested or first).of_px(px)
        return FigureLength(round(converted.value, _CONVERTED_DECIMALS), converted.unit).text

    return {"width": side(size.width, width_px), "height": side(size.height, height_px)}


def _export_root_size(drawing: Drawing) -> dict[str, str] | None:
    """The written root size of ``drawing``: ``export_root_size`` of the size it was bound to, if any."""

    size = getattr(drawing, _EXPORT_SIZE_ATTR, None)
    if not isinstance(size, RequestedFigureSize):
        return None
    return export_root_size(size, *drawing_size(drawing))


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
    """The text of one serialized drawing, held as its pieces (``DocumentParts``).

    ``prolog`` is what svgwrite writes before the root element of a file: the
    XML declaration and any stylesheet instructions. ``body`` is the root
    element alone; in-memory SVG, interactive SVG and the CairoSVG input use it
    without the prolog, as ``Drawing.tostring()`` gave them. It is joined the
    first time a reader asks for it: ``write`` writes the pieces in order, so a
    file that nothing else reads never holds the whole text (OV-3027).
    """

    parts: DocumentParts

    @property
    def prolog(self) -> str:
        return self.parts.prolog

    @cached_property
    def body(self) -> str:
        return self.parts.body()

    def write(self, handle: TextIO) -> None:
        """Write the file form, prolog then body, to an open text file; the same text as ``prolog + body``."""

        handle.write(self.prolog)
        for text in self.parts.texts():
            write_text(handle, text)


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

    def texts(self) -> Iterator[str]:
        """The root element's texts in order: head, each piece's texts, tail (the one order of the body)."""

        yield self.head
        yield from _piece_texts(self.pieces)
        yield self.tail

    def body(self) -> str:
        """The root element: ``texts`` joined. A reader that can take the texts in turn does not join them."""

        return "".join(self.texts())


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

    return SvgDocument(document_parts(drawing, root_attribs=_export_root_size(drawing)))


def raster_svg_body(drawing: Drawing, written: SvgDocument | None = None) -> str:
    """The text a PNG of ``drawing`` is drawn from: the root element with the drawing's own px size.

    ``written`` (``write_svg_document`` of ``drawing``) is that text when no requested size changed its root. A
    drawing fitted to a requested size writes its root in the requested unit, rounded (``export_root_size``), which
    can move a PNG side one pixel off ``png_pixel_size`` of the drawn px, the size the export check and the Web PNG
    take; its PNG is drawn from the px root instead. The viewBox and the content are the same text.
    """

    if written is not None and _export_root_size(drawing) is None:
        return written.body
    return document_parts(drawing).body()


__all__ = [
    "ContentPiece",
    "DocumentParts",
    "Piece",
    "RawPiece",
    "SvgDocument",
    "bind_export_size",
    "document_parts",
    "export_root_size",
    "raster_svg_body",
    "write_svg_document",
    "write_text",
]
