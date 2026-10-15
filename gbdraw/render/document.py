"""SVG document text: the one serialization of a drawing per request."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TextIO

from svgwrite import Drawing
from svgwrite.base import BaseElement

from gbdraw.render.group_content import ContentGroup, element_text, start_tag


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


_DECLARATION = '<?xml version="1.0" encoding="utf-8" ?>\n'
_STYLESHEET = (
    '<?xml-stylesheet href="%s" type="text/css" title="%s" alternate="%s" media="%s"?>\n'
)


def _child_text(element: BaseElement) -> str:
    if (
        isinstance(element, ContentGroup)
        and element.content.children
        and element.holds_content_children()
    ):
        return start_tag(element) + element.content.children_text + f"</{element.elementname}>"
    # defs, a group without content, and an empty group (``<g ... />``).
    return element_text(element)


def write_svg_document(drawing: Drawing) -> SvgDocument:
    """Serialize ``drawing`` once, as ``Drawing.write`` would write it.

    The prolog is svgwrite's: the XML declaration and the stylesheet
    instructions. The body is the root start tag, then each top-level child,
    then the root end tag. A placed content group with children is written as
    its start tag, the stored text of its children and its end tag; every
    other child by svgwrite's own serialization.
    """

    prolog = _DECLARATION + "".join(
        _STYLESHEET % stylesheet for stylesheet in drawing._stylesheets
    )
    if not drawing.elements:
        return SvgDocument(prolog=prolog, body=drawing.tostring())
    body = "".join(
        (
            start_tag(drawing),
            *(_child_text(element) for element in drawing.elements),
            f"</{drawing.elementname}>",
        )
    )
    return SvgDocument(prolog=prolog, body=body)


__all__ = ["SvgDocument", "write_svg_document"]
