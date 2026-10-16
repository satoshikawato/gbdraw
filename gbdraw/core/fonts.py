"""Which font file gbdraw measures a CSS ``font-family`` list with.

This module is the one owner of font availability. A family resolves, in this
order, to a font file registered in this process (``register_font_file``), to
a bundled Liberation face (by name or alias), to the file fontconfig
(``fc-match``) finds for exactly that family, or, when fontconfig is missing or
does not answer, to a file of that family in the platform's font folders. The
first family of a list that resolves is the one text is measured with; when
none resolves, Liberation Sans.

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
import sys
from dataclasses import dataclass
from importlib import resources
from typing import Callable, Dict, Iterable, Iterator, Literal, Optional, Sequence, Union

from fontTools.ttLib import TTCollection, TTFont

from gbdraw.core.numeric import ascii_float
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
# A cold fontconfig cache can take seconds to build on first use.
_FONTCONFIG_TIMEOUT_SECONDS = 10


@dataclass(frozen=True)
class FontFace:
    """One face of a font file.

    ``index`` is the face in a collection (.ttc), ``instance`` the 1-based
    named instance of a variable font (0: its default). ``content_id``
    identifies the bytes of a registered font; it is empty for bundled and
    system files, which do not change while gbdraw runs.
    """

    path: str
    index: int = 0
    instance: int = 0
    content_id: str = ""
    style: str = "Regular"

    @property
    def measurement_name(self) -> str:
        """The name text measurement reports for this face.

        A bundled or system face: its file stem ("LiberationSans-Regular"),
        with the collection index and named instance when they are not 0.
        A registered face: ``user-<16 hex of its SHA-256>-<style>``, the same
        for the same bytes whatever the file is called.
        """
        index = f"-{self.index}" if self.index else ""
        if self.content_id:
            return f"user-{self.content_id}{index}-{self.style}"
        instance = f"-i{self.instance}" if self.instance else ""
        return f"{os.path.splitext(os.path.basename(self.path))[0]}{index}{instance}"


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
    do not resolve, in their order; ``unchecked`` is the part of them that
    could not be looked up because fontconfig did not answer.
    """

    requested: str
    effective: str
    missing: tuple[str, ...]
    unchecked: tuple[str, ...] = ()

    @property
    def substituted(self) -> bool:
        """Whether ``effective`` differs from ``requested``."""
        return bool(self.missing)

    @property
    def replaced(self) -> tuple[str, ...]:
        """The requested families before the first one a drawing uses: what
        is drawn differs from the requested first family only when this is
        not empty (OV-743). A missing family after one that resolves is
        dropped from ``effective`` without changing the text."""
        missing = set(self.missing)
        leading: list[str] = []
        for family in parse_font_family_list(self.requested):
            if family not in missing:
                break
            leading.append(family)
        return tuple(leading)

    @property
    def nearest(self) -> str:
        """The family text is drawn with: the first named family of
        ``effective`` (``effective`` itself when it names only generics)."""
        return next((name for name in parse_font_family_list(self.effective) if not _is_generic(name)), self.effective)

    @property
    def message(self) -> Optional[str]:
        """The warning for a substitution: None when the requested first
        family is drawn (``replaced`` is empty)."""
        replaced = self.replaced
        if not replaced:
            return None
        parts: list[str] = []
        unchecked = [name for name in replaced if name in self.unchecked]
        unavailable = [name for name in replaced if name not in self.unchecked]
        if unavailable:
            noun, verb = ("Font", "is") if len(unavailable) == 1 else ("Fonts", "are")
            parts.append(f"{noun} {_css_family_list(unavailable)} {verb} not available")
        if unchecked:
            noun = "Font" if len(unchecked) == 1 else "Fonts"
            parts.append(
                f"{noun} {_css_family_list(unchecked)} could not be checked because fontconfig did not answer"
            )
        return (
            f"{'; '.join(parts)}; text uses {css_family_token(self.nearest)} instead "
            f"(font-family: {self.effective})."
        )


# ------------------------------------------------------------------
#  Names, weights and styles
# ------------------------------------------------------------------
def _normalize_family_name(value: str) -> str:
    return value.strip().strip("\"'").strip()


def _comparable_name(value: str) -> str:
    return _normalize_family_name(value).lower().replace(" ", "")


def _is_generic(family: str) -> bool:
    return _normalize_family_name(family).lower() in _CSS_GENERIC_FAMILIES


def css_family_token(family: str) -> str:
    """One family as a CSS token: a generic keyword bare, a name quoted."""
    name = _normalize_family_name(family)
    if name.lower() in _CSS_GENERIC_FAMILIES:
        return name.lower()
    quote = '"' if "'" in name else "'"
    return f"{quote}{name}{quote}"


def _css_family_list(families: Sequence[str]) -> str:
    return ", ".join(css_family_token(name) for name in families)


def parse_font_family_list(font_family: str) -> list[str]:
    """Split a CSS-style ``font-family`` list into unquoted family names."""
    return [
        family
        for family in (_normalize_family_name(part) for part in str(font_family).split(","))
        if family
    ]


def _css_weight(font_weight: FontWeight) -> int:
    """A CSS ``font-weight`` as a number (normal 400, bold and bolder 700)."""
    value = "normal" if font_weight is None else str(font_weight).strip().lower()
    if value in {"bold", "bolder"}:
        return 700
    try:
        return int(ascii_float(value))
    except (TypeError, ValueError):
        return 400


def _is_italic(font_style: Optional[str]) -> bool:
    value = "normal" if font_style is None else str(font_style).strip().lower()
    return value in {"italic", "oblique"}


def _nearest_weight(available: Sequence[int], desired: int) -> int:
    """The weight CSS font matching picks (CSS Fonts 4, section 5.2)."""
    if desired in available:
        return desired
    lighter = sorted((weight for weight in available if weight < desired), reverse=True)
    heavier = sorted(weight for weight in available if weight > desired)
    if 400 <= desired <= 500:
        up_to_500 = [weight for weight in heavier if weight <= 500]
        order = up_to_500 + lighter + [weight for weight in heavier if weight > 500]
    elif desired < 400:
        order = lighter + heavier
    else:
        order = heavier + lighter
    return order[0]


_FaceCandidate = tuple[int, bool, FontFace]


def _match_face(candidates: Sequence[_FaceCandidate], weight: int, italic: bool) -> Optional[FontFace]:
    """The face of a family CSS font matching picks: the requested style
    first, then the nearest weight."""
    if not candidates:
        return None
    pool = [candidate for candidate in candidates if candidate[1] == italic] or list(candidates)
    chosen = _nearest_weight([candidate[0] for candidate in pool], weight)
    return next(candidate[2] for candidate in pool if candidate[0] == chosen)


def face_weight_and_italic(font: TTFont) -> tuple[int, bool]:
    """A face's CSS weight (OS/2 ``usWeightClass``) and whether it is italic."""
    os2 = font["OS/2"] if "OS/2" in font else None
    mac_style = int(font["head"].macStyle) if "head" in font else 0
    weight = int(os2.usWeightClass) if os2 is not None else (700 if mac_style & 1 else 400)
    italic = bool((os2 is not None and os2.fsSelection & 1) or mac_style & 2)
    return weight, italic


# ------------------------------------------------------------------
#  Font files
# ------------------------------------------------------------------
_font_cache: Dict[FontFace, TTFont] = {}


def named_instance(font: TTFont, instance: int) -> TTFont:
    """A variable font at its 1-based named instance."""
    from fontTools.varLib import instancer

    coordinates = font["fvar"].instances[instance - 1].coordinates
    return instancer.instantiateVariableFont(font, dict(coordinates))


def load_font(face: FontFace) -> Optional[TTFont]:
    """Load one face, cached; None (with a warning) when it cannot be read.

    A registered face is loaded at registration; once it is unregistered it
    is not read again from its path, which may hold other bytes by then.
    """
    if face in _font_cache:
        return _font_cache[face]
    if face.content_id:
        raise ValueError(f"Font {face.measurement_name} ({face.path}) is no longer registered.")
    try:
        font = TTFont(face.path, fontNumber=face.index)
        if face.instance:
            font = named_instance(font, face.instance)
    except Exception as exc:
        logger.warning("Failed to load font file %s: %s", face.path, exc)
        return None
    _font_cache[face] = font
    return font


def _forget_face(face: FontFace) -> None:
    font = _font_cache.pop(face, None)
    if font is not None:
        font.close()
    _registered_data.pop(face, None)


def face_source_bytes(face: FontFace) -> bytes:
    """The bytes of the file a face is in (a whole collection for a .ttc face).

    A registered face: the bytes read at registration. Raises ``ValueError``
    for a registered face that is no longer registered and ``OSError`` for an
    unreadable file.
    """
    if face.content_id:
        data = _registered_data.get(face)
        if data is None:
            raise ValueError(f"Font {face.measurement_name} ({face.path}) is no longer registered.")
        return data
    with open(face.path, "rb") as handle:
        return handle.read()


def is_bundled_face(face: FontFace) -> bool:
    """Whether a face is one of gbdraw's own (Liberation) files."""
    return not face.content_id and face.path in bundled_font_paths()


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


def _bundled_variant(prefix: str, weight: int, italic: bool) -> Optional[FontFace]:
    # Liberation has weights 400 and 700; CSS matching picks 700 from 600 up.
    bold = weight >= 600
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
_registered: Dict[str, list[RegisteredFont]] = {}
# The file bytes of each registered face (the BytesIO its TTFont reads shares them).
_registered_data: Dict[FontFace, bytes] = {}


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


def _read_font_file(font_path: str) -> bytes:
    try:
        with open(font_path, "rb") as handle:
            return handle.read()
    except OSError as exc:
        raise ValidationError(
            f"Cannot read font file {font_path}: {exc}",
            diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
        ) from exc


def register_font_file(path: Union[str, os.PathLike[str]], *, index: int = 0) -> RegisteredFont:
    """Make a TrueType or OpenType file resolvable by its family name.

    The family is the typographic family name (name ID 16), else the family
    name (ID 1); the weight is OS/2 ``usWeightClass``; the face is italic when
    OS/2 or ``head`` says so. ``index`` selects a face of a collection (.ttc).
    A family keeps every registered face; text of a weight and style is
    measured with the face CSS font matching picks (the requested style, then
    the nearest weight). Registering a face with the same family, weight and
    style again replaces the earlier one. A registered family wins over the
    bundled alias of the same name (a registered "Times New Roman" is not
    Liberation Serif).
    """
    font_path = os.fspath(path)
    return _register_faces([(_read_font_file(font_path), index, font_path)])[0]


def register_font_file_batch(paths: Sequence[Union[str, os.PathLike[str]]]) -> tuple[RegisteredFont, ...]:
    """``register_font_file`` for several files, all or none: every file is
    read and checked before the registry changes, so a failure leaves the
    earlier registrations, and the faces a later file would replace, as they
    were."""
    font_paths = [os.fspath(path) for path in paths]
    return _register_faces([(_read_font_file(font_path), 0, font_path) for font_path in font_paths])


def register_font_file_faces(path: Union[str, os.PathLike[str]]) -> tuple[RegisteredFont, ...]:
    """``register_font_file`` for every face of a file: each face of a
    collection (.ttc, .otc), else the one face."""
    font_path = os.fspath(path)
    data = _read_font_file(font_path)
    count = 1
    if data[:4] == b"ttcf":
        try:
            collection = TTCollection(io.BytesIO(data), lazy=True)
            count = len(collection.fonts)
            collection.close()
        except Exception as exc:
            raise ValidationError(
                f"Cannot read font file {font_path}: {exc}",
                diagnostic={"code": "INPUT_INVALID", "reason": "FONT_FILE"},
            ) from exc
    return _register_faces([(data, index, font_path) for index in range(count)])


@dataclass(frozen=True)
class _PreparedFace:
    """A face read and checked, not registered yet."""

    data: bytes
    font: TTFont
    registered: RegisteredFont


def _register_faces(sources: Sequence[tuple[bytes, int, str]]) -> tuple[RegisteredFont, ...]:
    """Register (bytes, face index, path) sources, all or none."""
    prepared: list[_PreparedFace] = []
    try:
        for data, index, font_path in sources:
            prepared.append(_prepare_face(data, index, font_path))
    except BaseException:
        for face in prepared:
            face.font.close()
        raise
    return tuple(_commit_face(face) for face in prepared)


def _prepare_face(data: bytes, index: int, font_path: str) -> _PreparedFace:
    # Loaded from the bytes read once, so a later change to the file cannot
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
        weight, italic = face_weight_and_italic(font)
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
    return _PreparedFace(data=data, font=font, registered=registered)


def _commit_face(prepared: _PreparedFace) -> RegisteredFont:
    registered, font, data = prepared.registered, prepared.font, prepared.data
    face, weight, italic = registered.face, registered.weight, registered.italic
    faces = _registered.setdefault(_comparable_name(registered.family), [])
    for replaced in [entry for entry in faces if (entry.weight, entry.italic) == (weight, italic)]:
        faces.remove(replaced)
        if replaced.face != face:
            _forget_face(replaced.face)
    if face in _font_cache:
        font.close()
    else:
        _font_cache[face] = font
    _registered_data[face] = data
    faces.append(registered)
    _clear_resolution_caches()
    return registered


_fontconfig_files: set[str] = set()


def show_font_file_to_cairo(face: FontFace, family: str) -> None:
    """Add a registered face's file to fontconfig's application fonts, once
    per file, so CairoSVG draws text with it even when the font is not
    installed. The CairoSVG export calls it for the registered faces a
    conversion draws with (gbdraw.render.font_embedding). Without
    fontconfig, those outputs draw an installed font of that name or a
    substitute, and a warning (once per file) says so. The browser draws no
    raster through cairo and is skipped."""
    font_path = face.path
    if sys.platform == "emscripten" or font_path in _fontconfig_files:
        return
    _fontconfig_files.add(font_path)
    import ctypes
    import ctypes.util

    added = 0
    # The soname first: it is the library cairo already loaded, and finding
    # it needs no search (``find_library`` runs a subprocess).
    for soname in ("libfontconfig.so.1", None):
        library_path = soname or ctypes.util.find_library("fontconfig")
        if not library_path:
            continue
        try:
            add_file = ctypes.CDLL(library_path).FcConfigAppFontAddFile
        except (OSError, AttributeError):
            continue
        add_file.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        add_file.restype = ctypes.c_int
        added = add_file(None, os.fsencode(os.path.abspath(font_path)))
        break
    if not added:
        logger.warning(
            "Font file %s could not be made available to CairoSVG (fontconfig is missing or refused it): "
            "PNG, PDF, EPS and PS output draw %s with an installed font of that name or a substitute.",
            font_path,
            css_family_token(family),
        )


def unregister_font_file(path: Union[str, os.PathLike[str]]) -> tuple[RegisteredFont, ...]:
    """Forget every registration of ``path``; returns what was removed."""
    font_path = os.fspath(path)
    removed: list[RegisteredFont] = []
    for family in list(_registered):
        faces = _registered[family]
        for registered in [entry for entry in faces if entry.face.path == font_path]:
            faces.remove(registered)
            removed.append(registered)
            _forget_face(registered.face)
        if not faces:
            del _registered[family]
    if removed:
        _clear_resolution_caches()
    return tuple(removed)


def clear_registered_fonts() -> None:
    """Forget every registered font file."""
    for faces in _registered.values():
        for registered in faces:
            _forget_face(registered.face)
    _registered.clear()
    _clear_resolution_caches()


def _clear_resolution_caches() -> None:
    resolve_named_font_face.cache_clear()
    resolve_font_face.cache_clear()


def registered_fonts() -> tuple[RegisteredFont, ...]:
    return tuple(registered for faces in _registered.values() for registered in faces)


def _registered_face(family: str, weight: int, italic: bool) -> Optional[FontFace]:
    faces = _registered.get(_comparable_name(family), [])
    return _match_face([(entry.weight, entry.italic, entry.face) for entry in faces], weight, italic)


# ------------------------------------------------------------------
#  Fontconfig
# ------------------------------------------------------------------
_FontconfigStatus = Literal["answered", "absent", "failed"]
# Answers only: a lookup that failed is tried again.
_fontconfig_cache: Dict[tuple[str, bool, bool], Optional[FontFace]] = {}


def _fontconfig_pattern(family: str, bold: bool, italic: bool) -> str:
    escaped = "".join(f"\\{char}" if char in "\\-:," else char for char in family)
    return f"{escaped}:weight={'bold' if bold else 'regular'}:slant={'italic' if italic else 'roman'}"


def _fontconfig_face(family: str, weight: int, italic: bool) -> tuple[Optional[FontFace], _FontconfigStatus]:
    """The file fontconfig finds for exactly ``family``; None for a fallback."""
    bold = weight >= 600
    key = (_comparable_name(family), bold, italic)
    if key in _fontconfig_cache:
        return _fontconfig_cache[key], "answered"
    fc_match = shutil.which("fc-match")
    if not fc_match:
        return None, "absent"
    try:
        # UTF-8 is fontconfig's output encoding on every platform; undecodable
        # path bytes survive as surrogates, so the path still opens.
        result = subprocess.run(
            [fc_match, "--format=%{family}\t%{file}\t%{index}", _fontconfig_pattern(family, bold, italic)],
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="surrogateescape",
            timeout=_FONTCONFIG_TIMEOUT_SECONDS,
        )
        fields = result.stdout.strip().split("\t")
        raw_index = int(fields[2]) if len(fields) > 2 and fields[2].strip() else 0
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        logger.debug("Fontconfig lookup failed for %r: %s", family, exc)
        return None, "failed"
    if result.returncode != 0:
        return None, "failed"
    face: Optional[FontFace] = None
    if len(fields) >= 2:
        names = {_comparable_name(name) for name in fields[0].split(",")}
        file_path = fields[1]
        if (
            _comparable_name(family) in names
            and file_path.lower().endswith(_FONT_FILE_SUFFIXES)
            and os.path.isfile(file_path)
        ):
            # fontconfig's index: the named instance of a variable font in the
            # high 16 bits, the face of a collection in the low 16 bits.
            face = FontFace(file_path, raw_index & 0xFFFF, raw_index >> 16)
    _fontconfig_cache[key] = face
    return face, "answered"


# ------------------------------------------------------------------
#  Platform font folders (no fontconfig)
# ------------------------------------------------------------------
_platform_index: Optional[Dict[str, list[_FaceCandidate]]] = None


def _platform_font_dirs() -> tuple[str, ...]:
    if sys.platform == "emscripten":
        return ()
    home = os.path.expanduser("~")
    if sys.platform == "win32":
        dirs = [os.path.join(os.environ.get("WINDIR", r"C:\Windows"), "Fonts")]
        local = os.environ.get("LOCALAPPDATA")
        if local:
            dirs.append(os.path.join(local, "Microsoft", "Windows", "Fonts"))
    elif sys.platform == "darwin":
        dirs = ["/System/Library/Fonts", "/Library/Fonts", os.path.join(home, "Library", "Fonts")]
    else:
        dirs = [
            "/usr/share/fonts",
            "/usr/local/share/fonts",
            os.path.join(home, ".local", "share", "fonts"),
            os.path.join(home, ".fonts"),
        ]
    return tuple(directory for directory in dirs if os.path.isdir(directory))


def _font_file_faces(path: str) -> Iterator[tuple[int, TTFont]]:
    if path.lower().endswith((".ttc", ".otc")):
        collection = TTCollection(path, lazy=True)
        try:
            yield from enumerate(collection.fonts)
        finally:
            collection.close()
    else:
        font = TTFont(path, lazy=True)
        try:
            yield 0, font
        finally:
            font.close()


def _family_names(font: TTFont) -> set[str]:
    names: set[str] = set()
    for record in font["name"].names:
        if record.nameID in (1, 16):
            try:
                names.add(_comparable_name(record.toUnicode()))
            except (UnicodeDecodeError, ValueError):
                continue
    return names


def _index_font_dirs(dirs: Sequence[str]) -> Dict[str, list[_FaceCandidate]]:
    """Every face in ``dirs`` by family name (name IDs 1 and 16, all languages)."""
    index: Dict[str, list[_FaceCandidate]] = {}
    for directory in dirs:
        for root, _subdirs, files in os.walk(directory):
            for file_name in sorted(files):
                if not file_name.lower().endswith(_FONT_FILE_SUFFIXES):
                    continue
                path = os.path.join(root, file_name)
                try:
                    for face_index, font in _font_file_faces(path):
                        weight, italic = face_weight_and_italic(font)
                        for name in _family_names(font):
                            index.setdefault(name, []).append((weight, italic, FontFace(path, face_index)))
                except Exception as exc:
                    logger.debug("Skipping unreadable font file %s: %s", path, exc)
    return index


def _platform_face(family: str, weight: int, italic: bool) -> Optional[FontFace]:
    global _platform_index
    if _platform_index is None:
        dirs = _platform_font_dirs()
        logger.info("Indexing the fonts in %s (once per run).", ", ".join(dirs) or "no font folders")
        _platform_index = _index_font_dirs(dirs)
    return _match_face(_platform_index.get(_comparable_name(family), []), weight, italic)


# ------------------------------------------------------------------
#  Resolution
# ------------------------------------------------------------------
_RegisteredLookup = Callable[[str, int, bool], Optional[FontFace]]


def _registered_or_bundled_face(
    family: str, weight: int, italic: bool, registered: _RegisteredLookup = _registered_face
) -> Optional[FontFace]:
    """A registered face of ``family``, else its bundled alias face (review
    R2: a registered family wins over the alias of the same name; a generic
    keyword is never a registered family). ``registered`` stands in for the
    registry when the Web parity vectors are generated."""
    if not _is_generic(family):
        face = registered(family, weight, italic)
        if face is not None:
            return face
    prefix = _BUNDLED_FONT_FAMILY_ALIASES.get(_comparable_name(family))
    return _bundled_variant(prefix, weight, italic) if prefix is not None else None


def _lookup_family(family: str, weight: int, italic: bool) -> tuple[Optional[FontFace], bool]:
    """The face of one family, and whether the lookup was left unchecked
    because fontconfig did not answer."""
    face = _registered_or_bundled_face(family, weight, italic)
    if face is not None or _is_generic(family) or _comparable_name(family) in _BUNDLED_FONT_FAMILY_ALIASES:
        return face, False
    face, status = _fontconfig_face(family, weight, italic)
    if face is None and status != "answered":
        # fontconfig covers the platform folders when it answers.
        face = _platform_face(family, weight, italic)
    return face, face is None and status == "failed"


def _resolve_family(family: str, weight: int, italic: bool) -> Optional[FontFace]:
    return _lookup_family(family, weight, italic)[0]


def _first_resolving_family(
    font_family: str, lookup: Callable[[str], Optional[FontFace]]
) -> tuple[Optional[str], Optional[FontFace]]:
    """The first family of a list that ``lookup`` resolves, and its face."""
    for family in parse_font_family_list(font_family):
        face = lookup(family)
        if face is not None:
            return family, face
    return None, None


@functools.lru_cache(maxsize=1024)
def resolve_named_font_face(
    font_family: str,
    font_weight: FontWeight = "normal",
    font_style: Optional[str] = "normal",
) -> tuple[Optional[str], Optional[FontFace]]:
    """The first family of ``font_family`` that resolves for this weight and
    style, and its face; ``(None, None)`` when none does. The family is the
    name a CSS ``@font-face`` rule for the face declares."""
    weight, italic = _css_weight(font_weight), _is_italic(font_style)
    return _first_resolving_family(font_family, lambda family: _resolve_family(family, weight, italic))


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
    face = resolve_named_font_face(font_family, font_weight, font_style)[1]
    if face is None:
        face = _bundled_variant(_FALLBACK_PREFIX, _css_weight(font_weight), _is_italic(font_style))
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


def _resolve_font_family_with(
    font_family: str, lookup: Callable[[str], tuple[bool, bool]]
) -> FontFamilyResolution:
    """``resolve_font_family`` with ``lookup`` answering, for a named
    (non-generic) family, whether it resolves and whether that could not be
    checked. The Web runs a port of this rule
    (gbdraw/web/js/services/font-families.js) on the vectors of
    tests/fixtures/font_family_resolution_vectors.json.
    """
    requested = str(font_family)
    kept: list[str] = []
    missing: list[str] = []
    unchecked: list[str] = []
    for family in parse_font_family_list(requested):
        if _is_generic(family):
            kept.append(family)
            continue
        available, family_unchecked = lookup(family)
        if available:
            kept.append(family)
            continue
        missing.append(family)
        if family_unchecked:
            unchecked.append(family)
    if not missing:
        return FontFamilyResolution(requested=requested, effective=requested, missing=())
    generics = [family for family in kept if _is_generic(family)]
    if len(generics) == len(kept):
        generic = _nearest_generic(generics, tuple(missing))
        if not generics:
            kept.append(generic)
        position = next(index for index, family in enumerate(kept) if _is_generic(family))
        kept.insert(position, _LIBERATION_FAMILY_BY_GENERIC[generic])
    return FontFamilyResolution(
        requested=requested,
        effective=_css_family_list(kept),
        missing=tuple(missing),
        unchecked=tuple(unchecked),
    )


def custom_font_families_among(font_family: str, registered_families: Iterable[str]) -> tuple[str, ...]:
    """The named families of a list that gbdraw does not bundle, where the
    named families are registered (the one custom-font rule; the Web app
    passes the families loaded in its page).

    A Session stores only their names (D-02), so another computer may draw
    them with a substitute. A registered family named like a bundled alias
    is one of them: here it draws with its own file, elsewhere with Liberation.
    """
    registered = {_comparable_name(family) for family in registered_families}
    return tuple(
        family
        for family in parse_font_family_list(font_family)
        if not _is_generic(family)
        and (
            _comparable_name(family) not in _BUNDLED_FONT_FAMILY_ALIASES
            or _comparable_name(family) in registered
        )
    )


def custom_font_families(font_family: str) -> tuple[str, ...]:
    """``custom_font_families_among`` the families registered in this process."""
    return custom_font_families_among(font_family, _registered)


def resolve_font_family(font_family: str) -> FontFamilyResolution:
    """Decide which requested families a drawing can use, and the nearest
    substitute when none of the named ones is available.

    Generic keywords are kept. A named family is kept when it resolves (bundled,
    registered, fontconfig, or a platform font folder). When a family is
    dropped and no named family remains, the Liberation family of the list's
    generic (else of a guess from the missing names: mono, code, courier,
    consol -> Mono; serif without sans, times, georgia, garamond, mincho,
    ming, song, roman -> Serif; else Sans) is inserted before the first
    generic, and that generic is added when the list has none. A list whose
    families all resolve comes back unchanged.
    """

    def lookup(family: str) -> tuple[bool, bool]:
        face, unchecked = _lookup_family(family, 400, False)
        return face is not None, unchecked

    return _resolve_font_family_with(font_family, lookup)


def resolve_font_family_among(font_family: str, registered_families: Iterable[str]) -> FontFamilyResolution:
    """``resolve_font_family`` where only the bundled fonts and the named
    registered families are available, without fontconfig: what the Web app
    (Pyodide, fonts registered in the page) resolves.
    """
    registered = {_comparable_name(family) for family in registered_families}
    return _resolve_font_family_with(
        font_family,
        lambda family: (
            _comparable_name(family) in _BUNDLED_FONT_FAMILY_ALIASES or _comparable_name(family) in registered,
            False,
        ),
    )


def font_family_rules() -> dict[str, object]:
    """The data of the family rules above, for the Web port
    (tools/generate_legend_font_metrics.py writes it to
    gbdraw/web/js/utils/font-family-rules.generated.js).
    """
    return {
        "bundledFamilyAliases": dict(sorted(_BUNDLED_FONT_FAMILY_ALIASES.items())),
        "genericFamilies": sorted(_CSS_GENERIC_FAMILIES),
        "liberationFamilyByGeneric": dict(sorted(_LIBERATION_FAMILY_BY_GENERIC.items())),
        "monospaceNameHints": list(_MONOSPACE_NAME_HINTS),
        "serifNameHints": list(_SERIF_NAME_HINTS),
        "fallbackFace": _FALLBACK_PREFIX,
    }


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
    return f"{css_family_token(family)}, {css_family_token(generic)}"


def complete_font_family(font_family: str) -> str:
    """Append the generic family to a single family name that resolves.

    A list, a generic keyword, or a name that does not resolve comes back
    unchanged, so a later substitution (``resolve_font_family``) decides it.
    """
    families = parse_font_family_list(font_family)
    if len(families) != 1 or "," in str(font_family) or _is_generic(families[0]):
        return str(font_family)
    face = _resolve_family(families[0], 400, False)
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
    "custom_font_families",
    "custom_font_families_among",
    "font_family_list_for",
    "css_family_token",
    "face_source_bytes",
    "face_weight_and_italic",
    "font_family_rules",
    "font_generic_family",
    "is_bundled_face",
    "load_font",
    "named_instance",
    "parse_font_family_list",
    "register_font_file",
    "register_font_file_batch",
    "register_font_file_faces",
    "registered_fonts",
    "resolve_font_face",
    "resolve_font_family",
    "resolve_font_family_among",
    "resolve_named_font_face",
    "show_font_file_to_cairo",
    "unregister_font_file",
]
