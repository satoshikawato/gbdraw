"""The fonts an exported file carries (D-05, D-13): SVG ``@font-face``
subsets, TrueType outlines for the Web PDF, and the faces CairoSVG draws with.

This module is the one owner of font embedding. ``face_embedding`` is the one
rule of what an output kind (``svg``, ``document``, ``raster``) may do with a
face, by its OS/2 ``fsType``; every path below asks it.

``svg_font_embedding`` finds the faces an SVG draws its text with (each
element's ``font-family``, ``font-weight`` and ``font-style``, resolved by
``gbdraw.core.fonts`` as text was measured) and the characters each face
draws, and returns ``@font-face`` rules with a subset of every face gbdraw
does not bundle, as ``data:`` URLs in the font's own flavour (TrueType or
OpenType). A face the output may not carry is left out and named in
``refused_families``. ``embedded_fonts_block`` builds the
``<defs id="gbdraw-embedded-fonts">`` block of those rules: the CLI and Python
writers insert it with ``embed_font_css``, the Web export receives it from the
Worker and only inserts it, and ``remove_embedded_font_css`` takes it out
again (a Session never stores font data, D-02).

``pdf_fonts`` chooses the face the Web PDF export draws each text run with
(the face the text was measured with), as TrueType outlines of the characters
it draws (``truetype_bytes_for_pdf``; CFF and CFF2 outlines are converted with
``cu2qu``), and falls back to the nearest Liberation face, by name, when a
custom face cannot be embedded. ``fonts_for_cairosvg`` does the same for the
CLI and Python PDF, EPS and PS, and shows the registered faces a conversion
uses to fontconfig.
"""

from __future__ import annotations

import base64
import functools
import hashlib
import io
import logging
import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Collection, Dict, Literal, Optional, Sequence, Union
from xml.sax.saxutils import escape

from fontTools.ttLib import TTFont

from gbdraw.core import fonts
from gbdraw.exceptions import ValidationError

logger = logging.getLogger(__name__)

EMBEDDED_FONTS_ID = "gbdraw-embedded-fonts"

_SVG_NAMESPACE = "{http://www.w3.org/2000/svg}"
_TEXT_CONTAINERS = frozenset({"text", "tspan", "textPath", "a"})
_FONT_PROPERTIES = ("font-family", "font-weight", "font-style")
# OS/2 fsType: bit 1 restricted-license embedding; bits 2 and 3 (preview &
# print, editable) allow it even when bit 1 is set.
_FS_TYPE_RESTRICTED = 0x0002
_FS_TYPE_ALLOWED = 0x0004 | 0x0008
# Bit 8: the font may be embedded only whole, never as a subset (D-09).
_FS_TYPE_NO_SUBSETTING = 0x0100
# Bit 9: only the font's bitmaps may be embedded, not its outlines.
_FS_TYPE_BITMAP_ONLY = 0x0200
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

EmbeddingKind = Literal["svg", "document", "raster"]
FaceEmbedding = Literal["subset", "whole", "refused"]


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


def face_embedding(font: TTFont, kind: EmbeddingKind) -> FaceEmbedding:
    """How an output of ``kind`` may carry ``font``, by its OS/2 ``fsType``
    (D-13); the one rule every export asks.

    ``raster`` (PNG): a subset whatever the license; rasterizing is not
    embedding, and the font never leaves the page or process. ``svg``
    (SVG, Interactive SVG): refused when the license allows only
    restricted-license embedding (0x0002 without 0x0004 or 0x0008) or only
    bitmaps (0x0200); whole when it forbids subsetting (0x0100, D-09); else a
    subset. ``document`` (PDF, EPS, PS, whose writers always subset): refused
    in all three of those cases; else a subset.
    """
    if kind == "raster":
        return "subset"
    fs_type = int(font["OS/2"].fsType) if "OS/2" in font else 0
    restricted = fs_type & _FS_TYPE_RESTRICTED and not fs_type & _FS_TYPE_ALLOWED
    if restricted or fs_type & _FS_TYPE_BITMAP_ONLY:
        return "refused"
    if fs_type & _FS_TYPE_NO_SUBSETTING:
        return "refused" if kind == "document" else "whole"
    return "subset"


def _face_data_url(face: fonts.FontFace, characters: set[str], *, whole: bool) -> str:
    """The face as a ``data:`` URL: a subset of ``characters``, or the whole
    face when its license forbids subsetting."""
    return _face_data_url_for(face, frozenset() if whole else frozenset(characters), whole)


# Keyed like ``_truetype_bytes``: a registered face's identity is its bytes,
# so an export that draws the same characters again reuses the subset and its
# base64 text. Bounded, so faces unregistered since are let go.
@functools.lru_cache(maxsize=16)
def _face_data_url_for(face: fonts.FontFace, characters: frozenset[str], whole: bool) -> str:
    font = TTFont(io.BytesIO(fonts.face_source_bytes(face)), fontNumber=face.index, recalcTimestamp=False)
    if face.instance:
        font = fonts.named_instance(font, face.instance)
    if not whole:
        _subset(font, characters)
    font.flavor = None
    buffer = io.BytesIO()
    font.save(buffer)
    mime, css_format = ("font/otf", "opentype") if font.sfntVersion == "OTTO" else ("font/ttf", "truetype")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f'url(data:{mime};base64,{encoded}) format("{css_format}")'


def _subset(font: TTFont, characters: Collection[str]) -> None:
    from fontTools import subset

    options = subset.Options()
    options.hinting = False
    options.notdef_outline = True
    options.name_IDs = [1, 2, 4, 6, 16, 17]
    # Keep the OS/2 Unicode ranges as they are: recalculating them fails on
    # a font that sets a reserved bit (fontTools 4: "expected 0 <= int <=
    # 122"), and the ranges do not choose glyphs (OV-744).
    options.prune_unicode_ranges = False
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


def svg_font_embedding(svg_text: str, kind: Literal["svg", "raster"] = "svg") -> FontEmbedding:
    """The ``@font-face`` rules for the custom faces ``svg_text`` draws with.

    Each face gets one rule with its own weight and style and a subset of the
    characters drawn with it, or the whole face, as ``face_embedding``
    decides for ``kind`` (``svg`` for an SVG file, ``raster`` for the SVG the
    Web PNG export draws, which takes the whole face when subsetting
    fails). Bundled Liberation faces are never embedded.
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
        use = face_embedding(font, kind)
        if use == "refused":
            if family not in refused:
                refused.append(family)
            continue
        weight, italic = fonts.face_weight_and_italic(font)
        try:
            src = _face_data_url(face, characters, whole=use == "whole")
        except Exception:
            if kind != "raster" or use == "whole":
                raise
            # The raster SVG never leaves the page: draw with the whole face.
            src = _face_data_url(face, characters, whole=True)
        rules.append(
            "@font-face { "
            f"font-family: {fonts.css_family_token(family)}; "
            f"font-weight: {weight}; "
            f"font-style: {'italic' if italic else 'normal'}; "
            f"src: {src}; "
            "}"
        )
    if not rules and not refused:
        return _NO_EMBEDDING
    return FontEmbedding(css="\n".join(rules), refused_families=tuple(refused))


def embedded_fonts_block(embedding: FontEmbedding) -> str:
    """The ``<defs id="gbdraw-embedded-fonts">`` block that carries
    ``embedding.css``, the first child of an exported SVG's root; "" when
    there is no CSS. The one builder of the block: ``embed_font_css`` inserts
    it, and the Web export inserts the block the Worker returns."""
    if not embedding.css:
        return ""
    return f'<defs id="{EMBEDDED_FONTS_ID}"><style type="text/css">{escape(embedding.css)}</style></defs>'


def embed_font_css(svg_text: str, embedding: FontEmbedding) -> str:
    """``svg_text`` with ``embedded_fonts_block`` first in the root element;
    unchanged when there is no CSS."""
    block = embedded_fonts_block(embedding)
    if not block:
        return svg_text
    match = _ROOT_START_TAG.search(svg_text)
    if match is None:
        raise ValueError("The SVG has no <svg> root element to embed fonts in.")
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


# Keyed by face and characters: a subset is small, and a registered face's
# identity is its bytes. Bounded, so faces unregistered since are let go.
@functools.lru_cache(maxsize=16)
def _truetype_bytes(face: fonts.FontFace, characters: Optional[frozenset[str]]) -> bytes:
    try:
        data = fonts.face_source_bytes(face)
        font = TTFont(io.BytesIO(data), fontNumber=face.index, recalcTimestamp=False)
        named = _ensure_postscript_name(font)
        if (
            characters is None and not face.index and not face.instance
            and data[:4] == b"\x00\x01\x00\x00" and named
        ):
            return data
        if "fvar" in font:
            from fontTools.varLib import instancer

            if face.instance:
                font = fonts.named_instance(font, face.instance)
            else:
                font = instancer.instantiateVariableFont(
                    font, {axis.axisTag: axis.defaultValue for axis in font["fvar"].axes}
                )
        if characters is not None:
            # Before cu2qu, which then converts only the glyphs drawn (S3).
            _subset(font, set(characters))
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


def _ensure_postscript_name(font: TTFont) -> bool:
    """Give ``font`` a PostScript name (name ID 6) when it has none, from its
    family and style; jsPDF reads it from every font a PDF embeds. Returns
    whether the font already had one."""
    names = font["name"]
    if names.getDebugName(6):
        return True
    family = names.getDebugName(16) or names.getDebugName(1) or "Font"
    style = names.getDebugName(17) or names.getDebugName(2) or "Regular"
    # PostScript names are printable ASCII without spaces and []{}<>()/%.
    postscript = re.sub(r"[^!-~]|[\[\]{}<>()/%]", "", f"{family}-{style}")[:63]
    names.setName(postscript, 6, 3, 1, 0x409)
    return False


def truetype_bytes_for_pdf(
    face_or_path: Union[fonts.FontFace, str, os.PathLike[str]],
    characters: Optional[Collection[str]] = None,
) -> bytes:
    """One face as a TrueType-outline font file, for a PDF to embed.

    A path names the first face of a file. With ``characters``, the face is
    first subset to them (so a large CFF font converts only the glyphs a PDF
    draws); without, TrueType faces come back as they are (one face of a
    collection, a variable font at its named instance or default). CFF and
    CFF2 outlines are converted to quadratic TrueType outlines with the same
    advance widths. Cached by face and characters (a registered face's
    identity is its bytes). Raises ``ValidationError`` (``FONT_FILE``) when
    the font cannot be read or converted.
    """
    face = face_or_path if isinstance(face_or_path, fonts.FontFace) else fonts.FontFace(os.fspath(face_or_path))
    return _truetype_bytes(face, None if characters is None else frozenset(characters))


PdfSubstitutionReason = Literal["license", "font-file"]


@dataclass(frozen=True)
class PdfFontFace:
    """The font a PDF text run is drawn with, as jsPDF registers it: the
    ``family`` (the face's measurement name without its style) and the
    TrueType ``file`` (one per face and character set, so it identifies the
    bytes)."""

    family: str
    file: str


@dataclass(frozen=True)
class PdfFontSubstitution:
    """A custom ``family`` the PDF draws with ``substitute`` instead, and why:
    its license forbids embedding, or its file cannot become TrueType."""

    family: str
    substitute: str
    reason: PdfSubstitutionReason


@dataclass(frozen=True)
class PdfFonts:
    """``faces``: one per run, in order. ``fonts``: the TrueType bytes of
    each file the caller does not hold yet. ``substitutions``: each custom
    family drawn with a Liberation face instead, once."""

    faces: tuple[PdfFontFace, ...]
    fonts: Dict[str, bytes]
    substitutions: tuple[PdfFontSubstitution, ...]


PdfFontRun = tuple[str, bool, bool, str]  # font-family list, bold, italic, text drawn


def _drawable(text: str) -> frozenset[str]:
    return frozenset(char for char in text if char >= " ")


def _pdf_font_face(face: fonts.FontFace, characters: Optional[frozenset[str]] = None) -> PdfFontFace:
    """A bundled face is sent whole, under its own file name; a custom face
    as a subset whose file name also names the characters."""
    name = face.measurement_name
    family = name.rsplit("-", 1)[0]
    if characters is None:
        return PdfFontFace(family=family, file=f"{name}.ttf")
    digest = hashlib.sha256("".join(sorted(characters)).encode("utf-8", "surrogatepass")).hexdigest()[:12]
    return PdfFontFace(family=family, file=f"{name}-{digest}.ttf")


def _nearest_bundled(font_family: str, weight: int, style: str) -> tuple[str, fonts.FontFace]:
    """The Liberation family and face a list falls back to without its
    custom families (the nearest font, as at Session load)."""
    resolution = fonts.resolve_font_family_among(font_family, ())
    face = fonts.resolve_font_face(resolution.effective, weight, style)
    assert face is not None
    return resolution.nearest, face


def _refused_in_document(face: fonts.FontFace) -> bool:
    if fonts.is_bundled_face(face):
        return False
    font = fonts.load_font(face)
    return font is not None and face_embedding(font, "document") == "refused"


def pdf_fonts(runs: Sequence[PdfFontRun], *, known: Collection[str]) -> PdfFonts:
    """The faces a PDF draws its text runs with, and the TrueType bytes of
    those not in ``known``.

    A run is drawn with the face its text was measured with
    (``resolve_font_face``). A custom face is sent as a subset of the
    characters all runs draw with it. A custom face a document may not carry
    (``face_embedding``), or that cannot become TrueType outlines, is replaced
    with the nearest Liberation face and named in ``substitutions``.
    """
    resolved: list[tuple[str, int, str, Optional[str], fonts.FontFace]] = []
    characters: Dict[fonts.FontFace, set[str]] = {}
    for font_family, bold, italic, text in runs:
        weight, style = (700 if bold else 400), ("italic" if italic else "normal")
        family, face = fonts.resolve_named_font_face(font_family, weight, style)
        if face is None:
            face = fonts.resolve_font_face(font_family, weight, style)
        assert face is not None
        resolved.append((font_family, weight, style, family, face))
        characters.setdefault(face, set()).update(_drawable(text))
    faces: list[PdfFontFace] = []
    data: Dict[str, bytes] = {}
    substitutions: list[PdfFontSubstitution] = []
    usable: Dict[fonts.FontFace, Optional[PdfSubstitutionReason]] = {}

    def send(face: fonts.FontFace, subset: Optional[frozenset[str]]) -> PdfFontFace:
        chosen = _pdf_font_face(face, subset)
        if chosen.file not in known and chosen.file not in data:
            data[chosen.file] = truetype_bytes_for_pdf(face, subset)
        return chosen

    for font_family, weight, style, family, face in resolved:
        bundled = fonts.is_bundled_face(face)
        subset = None if bundled else frozenset(characters[face])
        if face not in usable:
            usable[face] = "license" if _refused_in_document(face) else None
            if usable[face] is None:
                try:
                    send(face, subset)
                except ValidationError:
                    usable[face] = "font-file"
        reason = usable[face]
        if reason is None:
            faces.append(send(face, subset))
            continue
        substitute, fallback = _nearest_bundled(font_family, weight, style)
        substitution = PdfFontSubstitution(family or font_family, substitute, reason)
        if substitution not in substitutions:
            substitutions.append(substitution)
        faces.append(send(fallback, None))
    return PdfFonts(faces=tuple(faces), fonts=data, substitutions=tuple(substitutions))


# ------------------------------------------------------------------
#  CairoSVG: PNG, PDF, EPS and PS of the CLI and Python API
# ------------------------------------------------------------------
_CAIROSVG_KINDS: Dict[str, EmbeddingKind] = {"png": "raster", "pdf": "document", "eps": "document", "ps": "document"}
_SVG_NAMESPACE_URI = "http://www.w3.org/2000/svg"


def _set_font_family(element: ET.Element, font_family: str) -> None:
    """Set ``font-family`` as an attribute, dropping a ``style`` declaration
    of it, which would win over the attribute."""
    style = element.get("style")
    if style is not None and "font-family" in style:
        kept = [part for part in style.split(";") if part.partition(":")[0].strip() != "font-family"]
        element.set("style", ";".join(kept))
    element.set("font-family", font_family)


def fonts_for_cairosvg(svg_source: bytes, fmt: str) -> bytes:
    """The SVG CairoSVG converts to ``fmt``, with the fonts it may draw (D-13).

    Each registered face the text draws with is shown to fontconfig
    (``fonts.show_font_file_to_cairo``), so PNG, PDF, EPS and PS draw it
    when it is not installed. In PDF, EPS and PS, which carry the fonts they
    draw with, the text of a custom face a document may not carry
    (``face_embedding``) is set in the nearest Liberation face instead, with
    a warning naming both, as the Web PDF does. The bytes come back
    unchanged when nothing is rewritten.
    """
    kind = _CAIROSVG_KINDS.get(str(fmt).lower())
    if kind is None:
        return svg_source
    try:
        root = ET.fromstring(svg_source)
    except ET.ParseError:
        return svg_source
    substitutions: list[tuple[str, str]] = []
    checked: Dict[fonts.FontFace, bool] = {}

    def refused(family: Optional[str], face: Optional[fonts.FontFace]) -> bool:
        if family is None or face is None or fonts.is_bundled_face(face):
            return False
        if face not in checked:
            checked[face] = kind == "document" and _refused_in_document(face)
            if not checked[face] and face.content_id:
                fonts.show_font_file_to_cairo(face, family)
        return checked[face]

    def settle(element: ET.Element, inherited: _FontKey, in_text: bool, parent_rewritten: bool) -> bool:
        name = _local_name(element.tag)
        key = _font_properties(element, inherited)
        in_text = name == "text" or (in_text and name in _TEXT_CONTAINERS)
        rewritten = False
        changed = False
        if in_text:
            family, face = fonts.resolve_named_font_face(*key)
            if refused(family, face):
                resolution = fonts.resolve_font_family_among(key[0], ())
                _set_font_family(element, resolution.effective)
                if (str(family), resolution.nearest) not in substitutions:
                    substitutions.append((str(family), resolution.nearest))
                rewritten = changed = True
            elif parent_rewritten:
                _set_font_family(element, key[0])
                changed = True
        for child in element:
            changed = settle(child, key, in_text, rewritten) or changed
        return changed

    if not settle(root, ("", "normal", "normal"), False, False):
        return svg_source
    label = str(fmt).upper()
    for family, substitute in substitutions:
        logger.warning(
            "Font %s is not embedded in the %s: its license (OS/2 fsType) does not allow embedding. "
            "Its text is drawn in %s.",
            fonts.css_family_token(family),
            label,
            fonts.css_family_token(substitute),
        )
    try:
        return ET.tostring(root, encoding="utf-8", xml_declaration=True, default_namespace=_SVG_NAMESPACE_URI)
    except ValueError:
        return ET.tostring(root, encoding="utf-8", xml_declaration=True)


__all__ = [
    "EMBEDDED_FONTS_ID",
    "EmbeddingKind",
    "FaceEmbedding",
    "FontEmbedding",
    "PdfFontFace",
    "PdfFontRun",
    "PdfFontSubstitution",
    "PdfFonts",
    "pdf_fonts",
    "embed_font_css",
    "embedded_fonts_block",
    "face_embedding",
    "fonts_for_cairosvg",
    "remove_embedded_font_css",
    "svg_font_embedding",
    "truetype_bytes_for_pdf",
]
