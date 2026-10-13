"""The fonts an exported file carries (D-05): SVG ``@font-face`` subsets and
TrueType outlines for PDF.

This module is the one owner of font embedding. ``svg_font_embedding`` finds
the faces an SVG draws its text with (each element's ``font-family``,
``font-weight`` and ``font-style``, resolved by ``gbdraw.core.fonts`` as text
was measured) and the characters each face draws, and returns ``@font-face``
rules with a subset of every face gbdraw does not bundle, as ``data:`` URLs in
the font's own flavour (TrueType or OpenType). A face whose OS/2 ``fsType``
allows only restricted-license embedding is left out and named in
``refused_families``. Two entry points share it: the Web export inserts the
CSS itself, and the CLI and Python writers call ``embed_font_css``, which puts
it in a ``<defs>`` of its own that ``remove_embedded_font_css`` takes out
again (a Session never stores font data, D-02).

``truetype_bytes_for_pdf`` gives a face's glyphs as TrueType outlines, which
a PDF embeds; CFF and CFF2 outlines are converted with ``cu2qu``.
"""

from __future__ import annotations

import base64
import functools
import io
import logging
import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Dict, Optional, Union
from xml.sax.saxutils import escape

from fontTools.ttLib import TTFont

from gbdraw.core import fonts
from gbdraw.exceptions import ValidationError

EMBEDDED_FONTS_ID = "gbdraw-embedded-fonts"

_SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"
_TEXT_CONTAINERS = frozenset({"text", "tspan", "textPath", "a"})
_FONT_PROPERTIES = ("font-family", "font-weight", "font-style")
# OS/2 fsType: bit 1 restricted-license embedding; bits 2 and 3 (preview &
# print, editable) allow it even when bit 1 is set.
_FS_TYPE_RESTRICTED = 0x0002
_FS_TYPE_ALLOWED = 0x0004 | 0x0008
# The root <svg ...> start tag (attribute values may hold ">").
_ROOT_START_TAG = re.compile(r"<svg\b(?:[^>\"']|\"[^\"]*\"|'[^']*')*>")
_EMBEDDED_DEFS = re.compile(
    rf'<defs id="{EMBEDDED_FONTS_ID}"><style type="text/css">.*?</style></defs>', re.DOTALL
)


@dataclass(frozen=True)
class FontEmbedding:
    """``@font-face`` rules for the custom faces an SVG uses ("" for none),
    and the families whose license forbids embedding them, by name."""

    css: str
    refused_families: tuple[str, ...]


_NO_EMBEDDING = FontEmbedding(css="", refused_families=())

_FontKey = tuple[str, str, str]  # font-family, font-weight, font-style


def _local_name(tag: object) -> str:
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def _font_properties(element: ET.Element, inherited: _FontKey) -> _FontKey:
    values: Dict[str, str] = dict(zip(_FONT_PROPERTIES, inherited))
    for declaration in element.get("style", "").split(";"):
        name, _sep, value = declaration.partition(":")
        if name.strip() in values and value.strip():
            values[name.strip()] = value.strip()
    for name in _FONT_PROPERTIES:
        attribute = element.get(name)
        if attribute is not None and attribute.strip():
            values[name] = attribute.strip()
    return (values["font-family"], values["font-weight"], values["font-style"])


def _drawn_characters(root: ET.Element) -> Dict[_FontKey, set[str]]:
    """The characters each computed font draws, from the SVG's text elements."""
    drawn: Dict[_FontKey, set[str]] = {}

    def add(key: _FontKey, value: Optional[str]) -> None:
        if value:
            drawn.setdefault(key, set()).update(char for char in value if char >= " ")

    def walk(element: ET.Element, inherited: _FontKey, in_text: bool) -> None:
        name = _local_name(element.tag)
        key = _font_properties(element, inherited)
        in_text = name == "text" or (in_text and name in _TEXT_CONTAINERS)
        if in_text:
            add(key, element.text)
        for child in element:
            walk(child, key, in_text)
            if in_text:
                add(key, child.tail)

    walk(root, ("", "normal", "normal"), False)
    return drawn


def _subset_data_url(face: fonts.FontFace, characters: set[str]) -> str:
    from fontTools import subset

    font = TTFont(io.BytesIO(fonts.face_source_bytes(face)), fontNumber=face.index, recalcTimestamp=False)
    if face.instance:
        font = fonts.named_instance(font, face.instance)
    options = subset.Options()
    options.hinting = False
    options.notdef_outline = True
    options.name_IDs = [1, 2, 4, 6, 16, 17]
    subsetter = subset.Subsetter(options)
    subsetter.populate(unicodes={ord(char) for char in characters})
    # The subsetter reports each table it drops (for example a font editor's
    # private table) at WARNING; that is not a problem of the export.
    subset_log = logging.getLogger("fontTools.subset")
    level = subset_log.level
    subset_log.setLevel(logging.ERROR)
    try:
        subsetter.subset(font)
    finally:
        subset_log.setLevel(level)
    font.flavor = None
    buffer = io.BytesIO()
    font.save(buffer)
    mime, css_format = ("font/otf", "opentype") if font.sfntVersion == "OTTO" else ("font/ttf", "truetype")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f'url(data:{mime};base64,{encoded}) format("{css_format}")'


def svg_font_embedding(svg_text: str) -> FontEmbedding:
    """The ``@font-face`` rules for the custom faces ``svg_text`` draws with.

    Each face gets one rule with its own weight and style and a subset of the
    characters drawn with it. Bundled Liberation faces are never embedded.
    """
    drawn = _drawn_characters(ET.fromstring(svg_text))
    by_face: Dict[tuple[str, fonts.FontFace], set[str]] = {}
    for (font_family, font_weight, font_style), characters in drawn.items():
        family, face = fonts.resolve_named_font_face(font_family, font_weight, font_style)
        if family is None or face is None or fonts.is_bundled_face(face):
            continue
        by_face.setdefault((family, face), set()).update(characters)
    rules: list[str] = []
    refused: list[str] = []
    for (family, face), characters in sorted(by_face.items(), key=lambda item: (item[0][0], item[0][1].measurement_name)):
        font = fonts.load_font(face)
        if font is None:
            continue
        os2 = font["OS/2"] if "OS/2" in font else None
        fs_type = int(os2.fsType) if os2 is not None else 0
        if fs_type & _FS_TYPE_RESTRICTED and not fs_type & _FS_TYPE_ALLOWED:
            if family not in refused:
                refused.append(family)
            continue
        weight, italic = fonts.face_weight_and_italic(font)
        rules.append(
            "@font-face { "
            f"font-family: {fonts.css_family_token(family)}; "
            f"font-weight: {weight}; "
            f"font-style: {'italic' if italic else 'normal'}; "
            f"src: {_subset_data_url(face, characters)}; "
            "}"
        )
    if not rules and not refused:
        return _NO_EMBEDDING
    return FontEmbedding(css="\n".join(rules), refused_families=tuple(refused))


def embed_font_css(svg_text: str, embedding: FontEmbedding) -> str:
    """``svg_text`` with ``embedding.css`` in a ``<defs>`` of its own, first
    in the root element; unchanged when there is no CSS."""
    if not embedding.css:
        return svg_text
    match = _ROOT_START_TAG.search(svg_text)
    if match is None:
        raise ValueError("The SVG has no <svg> root element to embed fonts in.")
    block = f'<defs id="{EMBEDDED_FONTS_ID}"><style type="text/css">{escape(embedding.css)}</style></defs>'
    return svg_text[: match.end()] + block + svg_text[match.end():]


def remove_embedded_font_css(svg_text: str) -> str:
    """``svg_text`` without the block ``embed_font_css`` inserted."""
    return _EMBEDDED_DEFS.sub("", svg_text, count=1)


# ------------------------------------------------------------------
#  PDF: TrueType outlines
# ------------------------------------------------------------------
_POST_FORMAT = 2.0


def _to_quadratic(font: TTFont) -> None:
    """Replace CFF/CFF2 outlines with TrueType ``glyf`` outlines (the
    fontTools ``otf2ttf`` recipe); advance widths stay."""
    from fontTools.pens.cu2quPen import Cu2QuPen
    from fontTools.pens.ttGlyphPen import TTGlyphPen
    from fontTools.ttLib import newTable

    glyph_order = font.getGlyphOrder()
    glyph_set = font.getGlyphSet()
    quadratic = {}
    for name in glyph_order:
        pen = TTGlyphPen(glyph_set)
        glyph_set[name].draw(Cu2QuPen(pen, 1.0, reverse_direction=True))
        quadratic[name] = pen.glyph()
    font["loca"] = newTable("loca")
    font["glyf"] = glyf = newTable("glyf")
    glyf.glyphOrder = glyph_order
    glyf.glyphs = quadratic
    for tag in ("CFF ", "CFF2", "VORG"):
        if tag in font:
            del font[tag]
    glyf.compile(font)
    hmtx = font["hmtx"]
    for name, glyph in quadratic.items():
        if hasattr(glyph, "xMin"):
            hmtx[name] = (hmtx[name][0], glyph.xMin)
    font["maxp"] = maxp = newTable("maxp")
    maxp.tableVersion = 0x00010000
    maxp.maxZones = 1
    for field in (
        "maxTwilightPoints", "maxStorage", "maxFunctionDefs", "maxInstructionDefs",
        "maxStackElements", "maxSizeOfInstructions",
    ):
        setattr(maxp, field, 0)
    maxp.maxComponentElements = max(
        (len(glyph.components) for glyph in quadratic.values() if hasattr(glyph, "components")), default=0
    )
    maxp.compile(font)
    post = font["post"]
    post.formatType = _POST_FORMAT
    post.extraNames = []
    post.mapping = {}
    post.glyphOrder = glyph_order
    try:
        post.compile(font)
    except OverflowError:
        post.formatType = 3
    font.sfntVersion = "\x00\x01\x00\x00"


@functools.lru_cache(maxsize=16)
def _truetype_bytes(face: fonts.FontFace) -> bytes:
    try:
        data = fonts.face_source_bytes(face)
        if not face.index and not face.instance and data[:4] == b"\x00\x01\x00\x00":
            return data
        font = TTFont(io.BytesIO(data), fontNumber=face.index, recalcTimestamp=False)
        if "fvar" in font:
            from fontTools.varLib import instancer

            if face.instance:
                font = fonts.named_instance(font, face.instance)
            else:
                font = instancer.instantiateVariableFont(
                    font, {axis.axisTag: axis.defaultValue for axis in font["fvar"].axes}
                )
        if "glyf" not in font:
            if "CFF " not in font and "CFF2" not in font:
                raise ValueError("the font has neither TrueType nor CFF outlines")
            _to_quadratic(font)
        buffer = io.BytesIO()
        font.save(buffer)
        return buffer.getvalue()
    except Exception as exc:
        raise ValidationError(
            f"Cannot prepare font {face.path} for PDF: {exc}",
            diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
        ) from exc


def truetype_bytes_for_pdf(face_or_path: Union[fonts.FontFace, str, os.PathLike[str]]) -> bytes:
    """One face as a TrueType-outline font file, for a PDF to embed.

    A path names the first face of a file. TrueType faces come back as they
    are (one face of a collection, a variable font at its named instance or
    default); CFF and CFF2 outlines are converted to quadratic TrueType
    outlines with the same advance widths. Cached by face (a registered face's
    identity is its bytes). Raises ``ValidationError`` (``FONT_FILE``) when
    the font cannot be read or converted.
    """
    face = face_or_path if isinstance(face_or_path, fonts.FontFace) else fonts.FontFace(os.fspath(face_or_path))
    return _truetype_bytes(face)


__all__ = [
    "EMBEDDED_FONTS_ID",
    "FontEmbedding",
    "embed_font_css",
    "remove_embedded_font_css",
    "svg_font_embedding",
    "truetype_bytes_for_pdf",
]
