"""Which font file gbdraw measures a CSS ``font-family`` list with.

This module is the one owner of font availability. A family resolves, in this
order, to a bundled Liberation face (by name or alias), to a font file
registered in this process (``register_font_file``), or to the file fontconfig
(``fc-match``) finds for exactly that family. The first family of a list that
resolves is the one text is measured with; when none resolves, Liberation Sans.

``resolve_font_family`` turns a requested list into the list a drawing uses:
the families that resolve, in their order, plus the generic family, with the
nearest Liberation family inserted when no named family is left. Text is then
measured and drawn with the same font.
"""

from __future__ import annotations

import functools
import hashlib
import io
import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from importlib import resources
from typing import Dict, Optional, Union

from fontTools.ttLib import TTFont

from gbdraw.exceptions import ValidationError

logger = logging.getLogger(__name__)

FontWeight = Union[str, int, float, None]

_BUNDLED_FONT_FAMILY_ALIASES = {
    "liberationsans": "LiberationSans",
    "arial": "LiberationSans",
    "helvetica": "LiberationSans",
    "nimbussansl": "LiberationSans",
    "sans-serif": "LiberationSans",
    "sansserif": "LiberationSans",
    "liberationserif": "LiberationSerif",
    "timesnewroman": "LiberationSerif",
    "times": "LiberationSerif",
    "serif": "LiberationSerif",
    "liberationmono": "LiberationMono",
    "couriernew": "LiberationMono",
    "courier": "LiberationMono",
    "monospace": "LiberationMono",
}
_FALLBACK_PREFIX = "LiberationSans"
# CSS generic family keywords: kept in a list, never looked up as a font.
_CSS_GENERIC_FAMILIES = frozenset({
    "serif", "sans-serif", "monospace", "cursive", "fantasy", "system-ui", "ui-serif",
    "ui-sans-serif", "ui-monospace", "ui-rounded", "math", "emoji", "fangsong",
})
_LIBERATION_FAMILY_BY_GENERIC = {
    "serif": "Liberation Serif",
    "sans-serif": "Liberation Sans",
    "monospace": "Liberation Mono",
}
_MONOSPACE_NAME_HINTS = ("mono", "code", "courier", "consol")
_SERIF_NAME_HINTS = ("times", "georgia", "garamond", "mincho", "ming", "song", "roman")
_FONT_FILE_SUFFIXES = (".ttf", ".otf", ".ttc", ".otc")
_FONTCONFIG_TIMEOUT_SECONDS = 2


@dataclass(frozen=True)
class FontFace:
    """One face of a font file: the file and its index in a collection.

    ``content_id`` identifies the bytes of a registered font; it is empty for
    bundled and fontconfig files, which do not change while gbdraw runs.
    """

    path: str
    index: int = 0
    content_id: str = ""
    style: str = "Regular"

    @property
    def measurement_name(self) -> str:
        """The name text measurement reports for this face.

        A bundled or fontconfig face: its file stem ("LiberationSans-Regular").
        A registered face: ``user-<16 hex of its SHA-256>-<style>``, the same
        for the same bytes whatever the file is called.
        """
        if not self.content_id:
            return os.path.splitext(os.path.basename(self.path))[0]
        index = f"-{self.index}" if self.index else ""
        return f"user-{self.content_id}{index}-{self.style}"


@dataclass(frozen=True)
class RegisteredFont:
    """A font file made resolvable by its family name (``register_font_file``)."""

    family: str
    weight: int
    italic: bool
    generic: str
    face: FontFace


@dataclass(frozen=True)
class FontFamilyResolution:
    """A requested ``font-family`` list and the list a drawing uses.

    ``effective`` names only families that resolve, so text is measured with
    the font it is drawn with. ``missing`` lists the requested families that
    do not resolve, in their order.
    """

    requested: str
    effective: str
    missing: tuple[str, ...]

    @property
    def substituted(self) -> bool:
        return bool(self.missing)

    @property
    def message(self) -> Optional[str]:
        """The warning for a substitution; None when nothing was replaced."""
        if not self.missing:
            return None
        missing = ", ".join(_css_family_token(name) for name in self.missing)
        noun = "Font" if len(self.missing) == 1 else "Fonts"
        verb = "is" if len(self.missing) == 1 else "are"
        replacement = next(
            (name for name in parse_font_family_list(self.effective) if not _is_generic(name)),
            self.effective,
        )
        return (
            f"{noun} {missing} {verb} not available; text uses "
            f"{_css_family_token(replacement)} instead (font-family: {self.effective})."
        )


# ------------------------------------------------------------------
#  Names
# ------------------------------------------------------------------
def _normalize_family_name(value: str) -> str:
    return value.strip().strip("\"'").strip()


def _comparable_name(value: str) -> str:
    return _normalize_family_name(value).lower().replace(" ", "")


def _is_generic(family: str) -> bool:
    return _normalize_family_name(family).lower() in _CSS_GENERIC_FAMILIES


def _css_family_token(family: str) -> str:
    name = _normalize_family_name(family)
    if name.lower() in _CSS_GENERIC_FAMILIES:
        return name.lower()
    quote = '"' if "'" in name else "'"
    return f"{quote}{name}{quote}"


def parse_font_family_list(font_family: str) -> list[str]:
    """Split a CSS-style ``font-family`` list into unquoted family names."""
    return [
        family
        for family in (_normalize_family_name(part) for part in str(font_family).split(","))
        if family
    ]


def _is_bold(font_weight: FontWeight) -> bool:
    value = "normal" if font_weight is None else str(font_weight).strip().lower()
    if value in {"bold", "bolder"}:
        return True
    try:
        return int(float(value)) >= 600
    except (TypeError, ValueError):
        return False


def _is_italic(font_style: Optional[str]) -> bool:
    value = "normal" if font_style is None else str(font_style).strip().lower()
    return value in {"italic", "oblique"}


# ------------------------------------------------------------------
#  Font files
# ------------------------------------------------------------------
_font_cache: Dict[FontFace, TTFont] = {}


def load_font(face: FontFace) -> Optional[TTFont]:
    """Load one face, cached; None (with a warning) when it cannot be read."""
    if face in _font_cache:
        return _font_cache[face]
    try:
        font = TTFont(face.path, fontNumber=face.index)
    except Exception as exc:
        logger.warning("Failed to load font file %s: %s", face.path, exc)
        return None
    _font_cache[face] = font
    return font


def _forget_face(face: FontFace) -> None:
    font = _font_cache.pop(face, None)
    if font is not None:
        font.close()


def _bundled_font_dirs():
    data_dir = resources.files("gbdraw.data")
    return (data_dir.joinpath("fonts"), data_dir)


def bundled_font_paths() -> tuple[str, ...]:
    """The bundled font files, in directory order."""
    paths: list[str] = []
    for font_dir in _bundled_font_dirs():
        try:
            if not font_dir.is_dir():
                continue
            for font_file in font_dir.iterdir():
                font_path = str(font_file)
                if font_file.name.lower().endswith((".ttf", ".otf")) and font_path not in paths:
                    paths.append(font_path)
        except Exception as exc:
            logger.debug("Bundled font directory lookup failed for %s: %s", font_dir, exc)
    return tuple(paths)


_bundled_by_stem: Optional[Dict[str, FontFace]] = None


def _bundled_faces() -> Dict[str, FontFace]:
    global _bundled_by_stem
    if _bundled_by_stem is None:
        _bundled_by_stem = {}
        for path in bundled_font_paths():
            face = FontFace(path)
            _bundled_by_stem.setdefault(face.measurement_name.lower(), face)
    return _bundled_by_stem


def _bundled_variant(prefix: str, bold: bool, italic: bool) -> Optional[FontFace]:
    suffix = {
        (True, True): "BoldItalic",
        (True, False): "Bold",
        (False, True): "Italic",
        (False, False): "Regular",
    }[(bold, italic)]
    faces = _bundled_faces()
    return faces.get(f"{prefix}-{suffix}".lower()) or faces.get(f"{prefix}-Regular".lower())


# ------------------------------------------------------------------
#  Registered font files
# ------------------------------------------------------------------
_registered: Dict[str, Dict[tuple[bool, bool], RegisteredFont]] = {}


def _generic_of_font(font: TTFont) -> str:
    post = font["post"] if "post" in font else None
    if post is not None and post.isFixedPitch:
        return "monospace"
    panose = getattr(font["OS/2"], "panose", None) if "OS/2" in font else None
    if panose is not None and panose.bFamilyType == 2:
        if 11 <= panose.bSerifStyle <= 13:
            return "sans-serif"
        if 2 <= panose.bSerifStyle <= 10:
            return "serif"
    return "sans-serif"


def _read_font(source: Union[str, io.BytesIO], index: int, path: str) -> TTFont:
    try:
        return TTFont(source, fontNumber=index)
    except Exception as exc:
        raise ValidationError(
            f"Cannot read font file {path}: {exc}",
            diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
        ) from exc


def register_font_file(path: Union[str, os.PathLike[str]], *, index: int = 0) -> RegisteredFont:
    """Make a TrueType or OpenType file resolvable by its family name.

    The family is the typographic family name (name ID 16), else the family
    name (ID 1); the weight is OS/2 ``usWeightClass``; the face is italic when
    OS/2 or ``head`` says so. ``index`` selects a face of a collection (.ttc).
    Registering the same family, weight class and style again replaces the
    earlier file; a family without the requested bold or italic face is
    measured with its regular face.
    """
    font_path = os.fspath(path)
    try:
        with open(font_path, "rb") as handle:
            data = handle.read()
    except OSError as exc:
        raise ValidationError(
            f"Cannot read font file {font_path}: {exc}",
            diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
        ) from exc
    # Loaded from the bytes read here, so a later change to the file cannot
    # change what this registration measures.
    font = _read_font(io.BytesIO(data), index, font_path)
    try:
        name_table = font["name"] if "name" in font else None
        family = (
            (name_table.getDebugName(16) or name_table.getDebugName(1)) if name_table is not None else None
        )
        family = _normalize_family_name(family or "")
        if not family:
            raise ValidationError(
                f"Font file {font_path} has no family name.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
            )
        os2 = font["OS/2"] if "OS/2" in font else None
        mac_style = int(font["head"].macStyle) if "head" in font else 0
        weight = int(os2.usWeightClass) if os2 is not None else (700 if mac_style & 1 else 400)
        italic = bool((os2 is not None and os2.fsSelection & 1) or mac_style & 2)
        generic = _generic_of_font(font)
    except BaseException:
        font.close()
        raise
    bold = weight >= 600
    face = FontFace(
        font_path,
        index,
        content_id=hashlib.sha256(data).hexdigest()[:16],
        style={(True, True): "BoldItalic", (True, False): "Bold", (False, True): "Italic"}.get(
            (bold, italic), "Regular"
        ),
    )
    registered = RegisteredFont(family=family, weight=weight, italic=italic, generic=generic, face=face)
    slots = _registered.setdefault(_comparable_name(family), {})
    replaced = slots.get((bold, italic))
    if replaced is not None and replaced.face != face:
        _forget_face(replaced.face)
    if face in _font_cache:
        font.close()
    else:
        _font_cache[face] = font
    slots[(bold, italic)] = registered
    resolve_font_face.cache_clear()
    return registered


def unregister_font_file(path: Union[str, os.PathLike[str]]) -> tuple[RegisteredFont, ...]:
    """Forget every registration of ``path``; returns what was removed."""
    font_path = os.fspath(path)
    removed: list[RegisteredFont] = []
    for family in list(_registered):
        slots = _registered[family]
        for slot, registered in list(slots.items()):
            if registered.face.path == font_path:
                removed.append(slots.pop(slot))
                _forget_face(registered.face)
        if not slots:
            del _registered[family]
    if removed:
        resolve_font_face.cache_clear()
    return tuple(removed)


def clear_registered_fonts() -> None:
    """Forget every registered font file."""
    for slots in _registered.values():
        for registered in slots.values():
            _forget_face(registered.face)
    _registered.clear()
    resolve_font_face.cache_clear()


def registered_fonts() -> tuple[RegisteredFont, ...]:
    return tuple(registered for slots in _registered.values() for registered in slots.values())


def _registered_face(family: str, bold: bool, italic: bool) -> Optional[FontFace]:
    slots = _registered.get(_comparable_name(family))
    if not slots:
        return None
    for slot in ((bold, italic), (False, italic), (bold, False), (False, False)):
        if slot in slots:
            return slots[slot].face
    return next(iter(slots.values())).face


# ------------------------------------------------------------------
#  Fontconfig
# ------------------------------------------------------------------
_fontconfig_cache: Dict[tuple[str, bool, bool], Optional[FontFace]] = {}


def _fontconfig_pattern(family: str, bold: bool, italic: bool) -> str:
    escaped = "".join(f"\\{char}" if char in "\\-:," else char for char in family)
    return f"{escaped}:weight={'bold' if bold else 'regular'}:slant={'italic' if italic else 'roman'}"


def _fontconfig_face(family: str, bold: bool, italic: bool) -> Optional[FontFace]:
    """The file fontconfig finds for exactly ``family``; None for a fallback."""
    key = (_comparable_name(family), bold, italic)
    if key in _fontconfig_cache:
        return _fontconfig_cache[key]
    face: Optional[FontFace] = None
    fc_match = shutil.which("fc-match")
    if fc_match:
        try:
            result = subprocess.run(
                [fc_match, "--format=%{family}\t%{file}\t%{index}", _fontconfig_pattern(family, bold, italic)],
                check=False,
                capture_output=True,
                text=True,
                timeout=_FONTCONFIG_TIMEOUT_SECONDS,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            logger.debug("Fontconfig lookup failed for %r: %s", family, exc)
        else:
            fields = result.stdout.strip().split("\t")
            if result.returncode == 0 and len(fields) >= 2:
                names = {_comparable_name(name) for name in fields[0].split(",")}
                file_path = fields[1]
                index = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
                if (
                    _comparable_name(family) in names
                    and file_path.lower().endswith(_FONT_FILE_SUFFIXES)
                    and os.path.isfile(file_path)
                ):
                    face = FontFace(file_path, index)
    _fontconfig_cache[key] = face
    return face


# ------------------------------------------------------------------
#  Resolution
# ------------------------------------------------------------------
def _resolve_family(family: str, bold: bool, italic: bool) -> Optional[FontFace]:
    prefix = _BUNDLED_FONT_FAMILY_ALIASES.get(_comparable_name(family))
    if prefix is not None:
        return _bundled_variant(prefix, bold, italic)
    face = _registered_face(family, bold, italic)
    if face is not None or _is_generic(family):
        return face
    return _fontconfig_face(family, bold, italic)


@functools.lru_cache(maxsize=1024)
def resolve_font_face(
    font_family: str,
    font_weight: FontWeight = "normal",
    font_style: Optional[str] = "normal",
) -> Optional[FontFace]:
    """The face text in ``font_family`` is measured with.

    The first family of the list that resolves wins; when none does, the
    Liberation Sans face of the weight and style. None only when no bundled
    font exists.
    """
    bold, italic = _is_bold(font_weight), _is_italic(font_style)
    face: Optional[FontFace] = None
    for family in parse_font_family_list(font_family):
        face = _resolve_family(family, bold, italic)
        if face is not None:
            break
    if face is None:
        face = _bundled_variant(_FALLBACK_PREFIX, bold, italic)
    if face is None:
        bundled = bundled_font_paths()
        face = FontFace(bundled[0]) if bundled else None
    return face


def _nearest_generic(generics: list[str], missing: tuple[str, ...]) -> str:
    for generic in generics:
        if generic.lower() in _LIBERATION_FAMILY_BY_GENERIC:
            return generic.lower()
    names = " ".join(missing).lower()
    if any(hint in names for hint in _MONOSPACE_NAME_HINTS):
        return "monospace"
    if ("serif" in names and "sans" not in names) or any(hint in names for hint in _SERIF_NAME_HINTS):
        return "serif"
    return "sans-serif"


def resolve_font_family(font_family: str) -> FontFamilyResolution:
    """Decide which requested families a drawing can use, and the nearest
    substitute when none of the named ones is available.

    Generic keywords are kept. A named family is kept when it resolves (bundled,
    registered, or fontconfig). When a family is dropped and no named family
    remains, the Liberation family of the list's generic (else of a guess from
    the missing names: mono, code, courier, consol -> Mono; serif without sans,
    times, georgia, garamond, mincho, ming, song, roman -> Serif; else Sans) is
    inserted before the first generic, and that generic is added when the list
    has none. A list whose families all resolve comes back unchanged.
    """
    requested = str(font_family)
    kept: list[str] = []
    missing: list[str] = []
    for family in parse_font_family_list(requested):
        if _is_generic(family) or _resolve_family(family, False, False) is not None:
            kept.append(family)
        else:
            missing.append(family)
    if not missing:
        return FontFamilyResolution(requested=requested, effective=requested, missing=())
    generics = [family for family in kept if _is_generic(family)]
    if len(generics) == len(kept):
        generic = _nearest_generic(generics, tuple(missing))
        if not generics:
            kept.append(generic)
        position = next(index for index, family in enumerate(kept) if _is_generic(family))
        kept.insert(position, _LIBERATION_FAMILY_BY_GENERIC[generic])
    effective = ", ".join(_css_family_token(family) for family in kept)
    return FontFamilyResolution(requested=requested, effective=effective, missing=tuple(missing))


# ------------------------------------------------------------------
#  Generic family of a font file
# ------------------------------------------------------------------
def font_generic_family(path: Union[str, os.PathLike[str]], *, index: int = 0) -> str:
    """``serif``, ``sans-serif`` or ``monospace`` for a font file.

    Fixed pitch (``post.isFixedPitch``) is monospace; a Latin-text PANOSE
    (family type 2) with serif style 11-13 is sans-serif and 2-10 serif;
    anything else is sans-serif.
    """
    font_path = os.fspath(path)
    font = _read_font(font_path, index, font_path)
    try:
        return _generic_of_font(font)
    finally:
        font.close()


def font_family_list_for(family: str, generic: str) -> str:
    """The stored ``font-family`` list for one family: ``'Name', generic``."""
    return f"{_css_family_token(family)}, {_css_family_token(generic)}"


def complete_font_family(font_family: str) -> str:
    """Append the generic family to a single family name that resolves.

    A list, a generic keyword, or a name that does not resolve comes back
    unchanged, so a later substitution (``resolve_font_family``) decides it.
    """
    families = parse_font_family_list(font_family)
    if len(families) != 1 or "," in str(font_family) or _is_generic(families[0]):
        return str(font_family)
    face = _resolve_family(families[0], False, False)
    if face is None:
        return str(font_family)
    font = load_font(face)
    generic = _generic_of_font(font) if font is not None else "sans-serif"
    return font_family_list_for(families[0], generic)


__all__ = [
    "FontFace",
    "FontFamilyResolution",
    "RegisteredFont",
    "bundled_font_paths",
    "clear_registered_fonts",
    "complete_font_family",
    "font_family_list_for",
    "font_generic_family",
    "load_font",
    "parse_font_family_list",
    "register_font_file",
    "registered_fonts",
    "resolve_font_face",
    "resolve_font_family",
    "unregister_font_file",
]
