"""Compose a recorded walkthrough into a captioned, camera-guided MP4.

Inputs are the run-local files written by ``walkthrough.record_walkthrough``.
This module never starts a browser, so captions, pacing, and camera cues can be
revised from one recording.
"""

from __future__ import annotations

import bisect
import functools
import hashlib
import json
import math
import os
import re
import subprocess
from dataclasses import dataclass, replace
from html.parser import HTMLParser
from pathlib import Path

from fontTools.ttLib import TTCollection, TTFont
from PIL import Image, ImageChops, ImageDraw, ImageFilter, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[3]
INTER = REPO_ROOT / "gbdraw/web/vendor/fonts/inter"
LOGO = REPO_ROOT / "gbdraw/web/assets/gbdraw-logo-title.png"
TRANSLATIONS = REPO_ROOT / "docs/videos/meet-gbdraw"
GALLERY_EXAMPLES = REPO_ROOT / "gbdraw/web/gallery/examples.json"
NOTO_CJK = Path("/usr/share/fonts/opentype/noto")
FONT_CACHE = Path.home() / ".cache/gbdraw-video-fonts"

FPS = 30
RADIUS = 18
INTRO_SECONDS = 3.4
OPENING_INTRO_SECONDS = 2.0   # the title card after an opening of highlights
OUTRO_SECONDS = 5.0
FADE_SECONDS = 0.5
FINALE_TEXT = ("Vector stays sharp at any zoom", "The downloaded SVG, magnified 7×")
MONTAGE_CAPTION = (None, "Microbial and organelle genomes",
                   "Highlights from a live session in the gbdraw web app, sped up")
INTRO_LINES = ("From a GenBank file to a publication-ready genome map",
               "A live session in the gbdraw web app. Waits are fast-forwarded and marked.")
TUTORIAL_NOTE = "A Gallery tutorial, recorded live in the gbdraw web app. Waits are fast-forwarded and marked."
INK = (241, 245, 249)
MUTED = (148, 170, 196)
ACCENT = (59, 130, 246)
PORTRAIT_MAX_ZOOM = 2.6
PORTRAIT_PAD = 1.2  # portrait: margin around a framed target
PORTRAIT_ZOOM_SLACK = 1.35  # portrait: a zoomed tall target shows at most this much more height
PROGRESS_GAP = 22  # portrait: step bar distance below the window
INSET_MARGIN = 18  # overview inset: distance from the window's edges
INSET_FADE = 0.5  # seconds for the inset to fade out and in at its other corner
SLANT = 0.2  # synthetic oblique for <i> runs; the app vendors no italic face


@dataclass(frozen=True)
class Layout:
    """Where the app window, captions, and brand sit in one output format."""

    name: str
    size: tuple[int, int]
    card: tuple[int, int, int, int]      # x, y, width, height of the app window
    portrait: bool
    title_size: int
    subtitle_size: int
    ui_scale: float                      # toasts and badges
    finale: str                          # directory of the vector-zoom frames
    audio: bool                          # add a silent track for social uploads
    inset: int = 0                       # width of the whole-stage overview inset; 0 for none

    @property
    def aspect(self) -> float:
        return self.card[2] / self.card[3]

    @property
    def caption_steps(self) -> tuple[int, int]:
        """The line pitch of a caption's title and of its subtitle."""

        return round(self.title_size * 1.18), round(self.subtitle_size * 1.35)

    def caption_rows(self, title_lines: int, subtitle_lines: int) -> tuple[int, int, int]:
        """Where a caption's title and subtitle start in its layer, and the layer's height.

        Landscape captions sit in one fixed band below the window, so their
        height is the band's; a portrait caption is as tall as its lines.
        """

        title_step, subtitle_step = self.caption_steps
        if self.portrait:
            return 0, title_step * title_lines, title_step * title_lines + subtitle_step * subtitle_lines + 12
        top = 18 if subtitle_lines else 30
        return top, top + title_step * (title_lines - 1) + 56, CAPTION_BAND

    def inset_box(self, viewport: tuple[int, int], corner: str) -> tuple[int, int, int, int]:
        """The overview inset in window pixels (x, y, width, height), at the ``"bottom"`` or ``"top"`` left corner."""

        height = round(self.inset * viewport[1] / viewport[0])
        y = INSET_MARGIN if corner == "top" else self.card[3] - INSET_MARGIN - height
        return INSET_MARGIN, y, self.inset, height

    def toast_origin(self, viewport: tuple[int, int], size: tuple[int, int], corner: str | None,
                     shown: float) -> tuple[int, int]:
        """Where a toast layer sits: centered at the window's foot, raised above the inset in the bottom corner."""

        x, y, width, height = self.card
        lift = 0
        if corner == "bottom":
            lift = round((height - self.inset_box(viewport, corner)[1] - INSET_MARGIN) * shown)
        return x + (width - size[0]) // 2, y + height - size[1] - 32 - lift

    def badge_origin(self, size: tuple[int, int]) -> tuple[int, int]:
        """The fast-forward badge's top right corner, clear of the inset on the left."""

        return self.card[0] + self.card[2] - size[0] - 22, self.card[1] + 20

    def free_band(self, viewport: tuple[int, int], corner: str | None) -> tuple[float, float]:
        """The window's rows, as fractions of its height, that the inset leaves free above or below it."""

        if not self.inset or corner is None:
            return 0.0, 1.0
        _, y, _, height = self.inset_box(viewport, corner)
        if corner == "top":
            return (y + height + INSET_MARGIN) / self.card[3], 1.0
        return 0.0, (y - INSET_MARGIN) / self.card[3]


CAPTION_BAND = 144   # px below the landscape window for a caption
LANDSCAPE = Layout("landscape", (1920, 1080), (160, 34, 1600, 900), False, 40, 26, 1.0, "finale", False)
# 9:16 for TikTok and RedNote. The window avoids the top 400 px and bottom
# 430 px, where the apps draw their own navigation, captions, and buttons.
# A tall window shows part of the app, so an inset of the whole stage marks
# where that part is.
PORTRAIT = Layout("portrait", (1080, 1920), (72, 430, 936, 1060), True, 54, 32, 1.3, "finale-portrait", True,
                  inset=340)
INSET_CORNERS = ("bottom", "top")   # the default corner first


@dataclass(frozen=True)
class Clip:
    """A cut of the recording from ``before`` s around mark ``first`` to ``after`` s around mark ``last``.

    Offsets are on the walkthrough timeline. A clip with a ``label`` shows it
    as its caption instead of the step's.
    """

    first: str
    before: float
    last: str
    after: float
    label: str = ""


@dataclass(frozen=True)
class Story:
    """How one journey's recording becomes videos: its intro, its outputs, and its language.

    The app shows ``language`` for the whole recording, and the videos show
    their text only in it: a take in another language is another journey.
    """

    intro: tuple[str, str]
    outputs: tuple[tuple[str, Layout, str], ...]   # (edit, layout, file stem without the language suffix)
    language: str = "en"
    outro: tuple[str, ...] = ()   # the outro's title and lines; empty keeps ``Composer.OUTRO_ROWS``
    # Highlights the walkthrough edit plays before a shorter title card, each under its label.
    opening: tuple[Clip, ...] = ()


TUTORIAL_IDS = ("HmmtDNA_ATskew", "tobacco-chloroplast", "BGC0000708-BGC0000713",
                "Vnig_TUMSAT-TG-2018", "hepatoplasmataceae_collinear",
                "vibrio-harveyi-group-collinear")
VERTICAL = "-vertical"   # a tutorial's take for RedNote and other phone feeds
# The UI languages each tutorial is recorded in. A journey in a language other
# than English finds every control by the app's own text for it (``Recorder.t``).
TUTORIAL_LANGUAGES = {name: ("en",) for name in TUTORIAL_IDS} | {
    "BGC0000708-BGC0000713": ("en", "zh-Hans"),
    "vibrio-harveyi-group-collinear": ("en", "zh-Hans"),
}


def take_name(example_id: str, *, vertical: bool = False, language: str = "en") -> str:
    """The journey name of one take: ``<id>[-vertical][.<language>]``, English without a suffix."""

    return example_id + (VERTICAL if vertical else "") + ("" if language == "en" else f".{language}")


def _tutorial(example_id: str, *, vertical: bool = False, language: str = "en") -> Story:
    """A Gallery tutorial opens with the Gallery example's own title.

    Its landscape take becomes the landscape video, and its portrait take,
    recorded for the tall window, the vertical one.
    """

    examples = json.loads(GALLERY_EXAMPLES.read_text(encoding="utf-8"))
    title = next(row["title"] for row in examples if row["id"] == example_id)
    stem = f"gbdraw-tutorial-{example_id}"
    output = ("walkthrough", PORTRAIT, stem + VERTICAL) if vertical else ("walkthrough", LANDSCAPE, stem)
    return Story((title, TUTORIAL_NOTE), (output,), language)


# The 0.14.0 promotion: one tall take in Simplified Chinese, without a sound track
# (``video.promo`` records it).
PROMO_ID = "promo-zh-0140-vertical"
PROMO_STORY = Story(
    ("gbdraw 0.14.0: what is new", "Microbial and organelle genome maps · now with a Simplified Chinese interface"),
    (("walkthrough", replace(PORTRAIT, audio=False), "gbdraw-0.14.0-vertical"),), "zh-Hans",
    outro=("gbdraw.app", "Command line: pip install gbdraw"),
    # The marks ``video.promo`` sets at the end of each result's hold.
    opening=(Clip("hl-lang", -2.0, "hl-lang", -0.5, "Simplified Chinese interface"),
             Clip("hl-font", -1.4, "hl-font", 0.0, "Chinese text in the figure"),
             Clip("hl-rotate", -1.4, "hl-rotate", 0.0, "Start at any base"),
             Clip("hl-lane", -1.4, "hl-lane", 0.0, "Features on lanes you choose"),
             Clip("hl-align", -1.4, "hl-align", 0.0, "Similarity-group alignment")))

# Keys are journey names; ``walkthrough.journeys()`` records the same set.
STORIES = {
    "meet-gbdraw": Story(INTRO_LINES, (
        ("walkthrough", LANDSCAPE, "gbdraw-walkthrough"),
        ("highlights", LANDSCAPE, "gbdraw-highlights"),
        ("highlights", PORTRAIT, "gbdraw-highlights-vertical"),
    )),
    **{take_name(name, vertical=vertical, language=language): _tutorial(name, vertical=vertical, language=language)
       for name in TUTORIAL_IDS for vertical in (False, True) for language in TUTORIAL_LANGUAGES[name]},
    PROMO_ID: PROMO_STORY,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _ease(value: float) -> float:
    value = min(1.0, max(0.0, value))
    return value * value * (3 - 2 * value)


def _contains(outer: dict, inner: dict, slack: float = 8) -> bool:
    return (outer["x"] - slack <= inner["x"] and outer["y"] - slack <= inner["y"]
            and inner["x"] + inner["width"] <= outer["x"] + outer["width"] + slack
            and inner["y"] + inner["height"] <= outer["y"] + outer["height"] + slack)


def _fonts(work: Path) -> dict[int, tuple[Path, int]]:
    """Convert the app's vendored Inter WOFF2 files once per render."""

    fonts = {}
    for weight in (400, 500, 600, 700):
        target = work / f"inter-{weight}.ttf"
        if not target.exists():
            font = TTFont(INTER / f"inter-latin-{weight}-normal.woff2")
            font.flavor = None
            font.save(target)
        fonts[weight] = (target, 0)
    return fonts


def _cjk_fonts() -> dict[int, tuple[Path, int]]:
    """Find Noto Sans SC: a font directory, the user cache, or fonts-noto-cjk."""

    roots = [Path(os.environ["GBDRAW_VIDEO_FONT_DIR"])] if os.environ.get("GBDRAW_VIDEO_FONT_DIR") else []
    faces: dict[int, tuple[Path, int]] = {}
    for root in roots + [FONT_CACHE]:
        for weight, name in ((400, "Regular"), (500, "Medium"), (700, "Bold")):
            path = root / f"NotoSansSC-{name}.otf"
            if path.is_file():
                faces[weight] = (path, 0)
        if 400 in faces:
            break
    else:
        for weight, name in ((400, "Regular"), (500, "Medium"), (700, "Bold")):
            path = NOTO_CJK / f"NotoSansCJK-{name}.ttc"
            if path.is_file():
                collection = TTCollection(str(path), lazy=True)
                index = next(i for i, font in enumerate(collection.fonts)
                             if "SC" in font["name"].getDebugName(1))
                faces[weight] = (path, index)
    if 400 not in faces:
        raise FileNotFoundError(
            "Chinese captions need Noto Sans SC: install fonts-noto-cjk, or put "
            "NotoSansSC-Regular/Medium/Bold.otf in GBDRAW_VIDEO_FONT_DIR or ~/.cache/gbdraw-video-fonts")
    faces.setdefault(500, faces[400])
    faces.setdefault(700, faces[500])
    faces[600] = faces[700]
    return faces


class Type:
    def __init__(self, fonts: dict[int, tuple[Path, int]]) -> None:
        self._fonts = fonts
        self._cache: dict = {}
        path, index = fonts[400]
        self._cmap = set(TTFont(str(path), fontNumber=index, lazy=True).getBestCmap())
        self.name = path.name

    def font(self, weight: int, size: int) -> ImageFont.FreeTypeFont:
        key = (weight, size)
        if key not in self._cache:
            path, index = self._fonts[weight]
            self._cache[key] = ImageFont.truetype(str(path), size, index=index)
        return self._cache[key]

    def check(self, text: str) -> None:
        missing = {char for char in text if ord(char) not in self._cmap and char not in "\n"}
        if missing:
            raise ValueError(f"{self.name} lacks glyphs {sorted(missing)!r} in {text!r}")

    def label_font(self, text: str, size: int) -> ImageFont.FreeTypeFont:
        """The regular face for app text such as a menu option; text it cannot show takes a system face that can.

        The Language menu names each language in its own script (ไทย, 한국어),
        which the caption face may lack; fontconfig picks a face for them.
        """

        missing = sorted({ord(char) for char in text if ord(char) not in self._cmap})
        if not missing:
            return self.font(400, size)
        key = ("fallback", tuple(missing), size)
        if key not in self._cache:
            charset = " ".join(f"{code:x}" for code in missing)
            try:
                found = subprocess.run(["fc-match", "-f", "%{file}|%{index}", f"sans-serif:charset={charset}"],
                                       capture_output=True, text=True, check=True).stdout
                path, index = found.rsplit("|", 1)
                self._cache[key] = ImageFont.truetype(path, size, index=int(index))
            except (OSError, ValueError, subprocess.CalledProcessError):
                self._cache[key] = self.font(400, size)
        return self._cache[key]


class Language:
    """On-screen text for one language; English is the recorded source text."""

    def __init__(self, code: str) -> None:
        self.code = code
        self.strings: dict[str, str] = {}
        self.patterns: list[tuple[re.Pattern, str]] = []
        if code != "en":
            data = json.loads((TRANSLATIONS / f"captions.{code}.json").read_text(encoding="utf-8"))
            self.strings = data["strings"]
            self.patterns = [(re.compile(row["match"]), row["replace"]) for row in data["patterns"]]

    @property
    def suffix(self) -> str:
        return "" if self.code == "en" else f".{self.code}"

    def type(self, work: Path) -> Type:
        return Type(_fonts(work) if self.code == "en" else _cjk_fonts())

    def __call__(self, text: str) -> str:
        if self.code == "en" or not text:
            return text
        if text in self.strings:
            return self.strings[text]
        for pattern, template in self.patterns:
            match = pattern.match(text)
            if match:
                return template.format(**match.groupdict())
        raise KeyError(f"No {self.code} text for {text!r}; add it to captions.{self.code}.json")


# -- timeline ---------------------------------------------------------------

@dataclass
class TimeMap:
    """Piecewise-linear map from recording time to output time."""

    starts: list[float]
    factors: list[float]
    outs: list[float]

    @classmethod
    def build(cls, events: list[dict], begin: float, end: float) -> "TimeMap":
        segments = []
        open_segment = None
        for event in events:
            if event["kind"] == "speed_start":
                open_segment = event
            elif event["kind"] == "speed_end" and open_segment is not None:
                length = event["t"] - open_segment["t"]
                factor = open_segment.get("factor") or max(1.0, length / open_segment["target"])
                kind = "fast" if open_segment.get("factor") else "wait"
                if factor > 1.05:
                    segments.append((open_segment["t"], event["t"], factor, kind))
                open_segment = None
        starts, factors, outs = [begin], [1.0], [0.0]
        for start, stop, factor, _ in segments:
            for moment, rate in ((start, factor), (stop, 1.0)):
                outs.append(outs[-1] + (moment - starts[-1]) / factors[-1])
                starts.append(moment)
                factors.append(rate)
        mapping = cls(starts, factors, outs)
        mapping.duration = mapping.out(end)
        mapping.segments = segments
        return mapping

    def out(self, moment: float) -> float:
        index = max(0, bisect.bisect_right(self.starts, moment) - 1)
        return self.outs[index] + (moment - self.starts[index]) / self.factors[index]

    def source(self, output: float) -> float:
        index = max(0, bisect.bisect_right(self.outs, output) - 1)
        return self.starts[index] + (output - self.outs[index]) * self.factors[index]

    def factor(self, output: float) -> float:
        return self.factors[max(0, bisect.bisect_right(self.outs, output) - 1)]

    def waiting(self, output: float) -> bool:
        """Whether ``output`` falls in a compressed wait (for example Generate)."""

        moment = self.source(output)
        return any(kind == "wait" and start <= moment < stop for start, stop, _, kind in self.segments)


class Camera:
    """Ease between logged framing targets, clamped to the page.

    A state is ``(cx, cy, h)``: the visible page center and height in CSS px.
    Landscape output follows the logged zoom. Portrait output re-frames each
    logged target (or the diagram, for whole-page shots) for a tall window,
    centered in the rows the overview inset leaves free. When the page edge
    would push a target under the inset, the inset moves to the other corner
    for that shot. While the stage shows only a smaller tab (a ``stage``
    event, for example an exported file opened in a tall tab), the window
    already shows that whole page, so the inset is hidden (corner None).
    """

    def __init__(self, events: list[dict], timemap: TimeMap, viewport: tuple[int, int], layout: Layout) -> None:
        self.viewport = viewport
        self.layout = layout
        self.aspect = layout.aspect
        self.max_h = min(viewport[1], viewport[0] / self.aspect)
        self.cues = []
        # The page rect the camera stays inside, from each output time on.
        whole = (0.0, 0.0, float(viewport[0]), float(viewport[1]))
        self._stage_starts, self._stages = [-math.inf], [whole]
        for event in events:
            if event["kind"] == "stage":
                rect = event["rect"]
                self._stage_starts.append(timemap.out(event["t"]))
                self._stages.append((rect["x"], rect["y"], rect["width"], rect["height"]))
            elif event["kind"] == "camera":
                self.cues.append((timemap.out(event["t"]), *self._target(event, self._stages[-1]),
                                  event["duration"]))
        self.starts = [cue[0] for cue in self.cues]
        self._states = []
        state = self._clamp(viewport[0] / 2, viewport[1] / 2, self.max_h)
        for index in range(len(self.cues)):
            if index:
                state = self._at(index - 1, self.starts[index])
            self._states.append(state)

    def _target(self, event: dict, stage: tuple[float, float, float, float]
                ) -> tuple[tuple[float, float, float], str | None]:
        """The cue's framing inside ``stage`` and the inset corner it leaves the target clear of (None: no inset)."""

        if not self.layout.portrait:
            return self._clamp(event["cx"], event["cy"], self.viewport[1] / max(1.0, event["zoom"]), stage), None
        rect = event.get("portrait") or event["rect"]
        subject = event.get("subject")
        # An overview too wide for the tall window shows the diagram instead;
        # a shot of part of it, such as a feature and its popup, keeps its target.
        if (rect["width"] * PORTRAIT_PAD / self.aspect > self.max_h and subject and not event.get("portrait")
                and _contains(rect, subject)):
            rect = subject
        if not (self.layout.inset and stage[2:] == tuple(self.viewport)):
            return self._frame(event, rect, None, stage), None
        framings = [(self._frame(event, rect, corner, stage), corner) for corner in INSET_CORNERS]
        # A target too large to clear the inset anywhere keeps the default corner.
        return next((row for row in framings if not self._covered(*row, rect)), framings[0])

    def _frame(self, event: dict, rect: dict, corner: str | None,
               stage: tuple[float, float, float, float]) -> tuple[float, float, float]:
        """Center the padded ``rect`` in the window rows that the inset at ``corner`` leaves free."""

        top, bottom = self.layout.free_band(self.viewport, corner)
        height = rect["height"] * PORTRAIT_PAD / (bottom - top)
        if event["zoom"] > 1 and not event.get("portrait"):
            # A tall target, such as the whole settings column, keeps close to
            # the logged zoom instead of shrinking to fit its full height.
            height = min(height, self.viewport[1] / event["zoom"] * PORTRAIT_ZOOM_SLACK)
        height = max(height, rect["width"] * PORTRAIT_PAD / self.aspect)
        height = min(max(height, self.max_h / PORTRAIT_MAX_ZOOM), self.max_h)
        cy = rect["y"] + rect["height"] / 2 + height * (0.5 - (top + bottom) / 2)
        return self._clamp(rect["x"] + rect["width"] / 2, cy, height, stage)

    def _covered(self, state: tuple[float, float, float], corner: str, rect: dict) -> bool:
        """Whether the inset at ``corner`` hides part of ``rect`` in the view ``state``."""

        cx, cy, height = state
        scale = self.layout.card[3] / height          # window px per page px
        left, top = cx - height * self.aspect / 2, cy - height / 2
        x, y, width, inset_height = self.layout.inset_box(self.viewport, corner)
        return ((rect["x"] - left) * scale < x + width and (rect["x"] + rect["width"] - left) * scale > x
                and (rect["y"] - top) * scale < y + inset_height and (rect["y"] + rect["height"] - top) * scale > y)

    def _clamp(self, cx: float, cy: float, height: float,
               stage: tuple[float, float, float, float] | None = None) -> tuple[float, float, float]:
        """Keep the view inside ``stage`` (x, y, width, height; default: the whole page)."""

        x, y, width, stage_h = stage or (0.0, 0.0, *self.viewport)
        height = min(height, self.max_h, stage_h, width / self.aspect)
        half_w, half_h = height * self.aspect / 2, height / 2
        return (min(max(cx, x + half_w), x + width - half_w),
                min(max(cy, y + half_h), y + stage_h - half_h), height)

    def _stage_at(self, output: float) -> tuple[float, float, float, float]:
        return self._stages[bisect.bisect_right(self._stage_starts, output) - 1]

    def _at(self, index: int, output: float) -> tuple[float, float, float]:
        start, target, _, duration = self.cues[index]
        origin = self._states[index]
        p = 1.0 if duration <= 0 else _ease((output - start) / duration)
        # Interpolate the visible height geometrically so zooming feels even.
        height = math.exp(math.log(origin[2]) + (math.log(target[2]) - math.log(origin[2])) * p)
        return (origin[0] + (target[0] - origin[0]) * p, origin[1] + (target[1] - origin[1]) * p, height)

    def at(self, output: float) -> tuple[float, float, float]:
        index = bisect.bisect_right(self.starts, output) - 1
        if index < 0:
            return self._clamp(self.viewport[0] / 2, self.viewport[1] / 2, self.max_h)
        return self._clamp(*self._at(index, output), self._stage_at(output))

    def inset_at(self, output: float) -> tuple[str | None, float]:
        """The inset's corner at ``output`` and its opacity: it fades out and back in when the corner changes."""

        index = bisect.bisect_right(self.starts, output) - 1
        if index < 0:
            return (INSET_CORNERS[0], 1.0) if self.layout.inset else (None, 0.0)
        start, _, corner, duration = self.cues[index]
        previous = self.cues[index - 1][2] if index else corner
        p = (output - start) / max(duration, INSET_FADE)
        if previous == corner or p >= 1:
            return corner, 1.0 if corner else 0.0
        if corner is None or previous is None:
            return (previous, 1 - p) if corner is None else (corner, p)
        return (previous, 1 - 2 * p) if p < 0.5 else (corner, 2 * p - 1)


class Pointer:
    def __init__(self, events: list[dict], timemap: TimeMap) -> None:
        moves = [(timemap.out(e["t"]), e["x"], e["y"]) for e in events if e["kind"] in ("move", "down", "up")]
        self.times = [row[0] for row in moves]
        self.points = [(row[1], row[2]) for row in moves]
        self.presses = [(timemap.out(e["t"]), e["x"], e["y"]) for e in events if e["kind"] == "down"]
        self.releases = [timemap.out(e["t"]) for e in events if e["kind"] == "up"]

    def at(self, output: float) -> tuple[float, float]:
        index = bisect.bisect_right(self.times, output)
        if index == 0:
            return self.points[0]
        if index >= len(self.times):
            return self.points[-1]
        t0, t1 = self.times[index - 1], self.times[index]
        (x0, y0), (x1, y1) = self.points[index - 1], self.points[index]
        p = 0 if t1 <= t0 else (output - t0) / (t1 - t0)
        return x0 + (x1 - x0) * p, y0 + (y1 - y0) * p

    def pressed(self, output: float) -> bool:
        index = bisect.bisect_right([p[0] for p in self.presses], output) - 1
        return index >= 0 and (index >= len(self.releases) or output < self.releases[index])

    def ripples(self, output: float) -> list[tuple[float, float, float]]:
        return [(x, y, (output - start) / 0.45) for start, x, y in self.presses if 0 <= output - start < 0.45]


# -- drawing ----------------------------------------------------------------

def _background(size: tuple[int, int]) -> Image.Image:
    width, height = size
    top, bottom = (15, 23, 42), (23, 43, 77)
    column = Image.new("RGB", (1, height))
    for y in range(height):
        p = y / (height - 1)
        column.putpixel((0, y), tuple(round(a + (b - a) * p) for a, b in zip(top, bottom)))
    image = column.resize(size)
    glow = Image.new("L", size, 0)
    ImageDraw.Draw(glow).ellipse((-width * 0.21, -height * 0.65, width * 0.73, height * 0.65), fill=70)
    glow = glow.filter(ImageFilter.GaussianBlur(220))
    return Image.composite(Image.new("RGB", size, (37, 99, 235)), image, glow.point(lambda v: v // 3))


def _card_layers(background: Image.Image, card: tuple[int, int, int, int]) -> tuple[Image.Image, Image.Image]:
    x, y, width, height = card
    shadow = Image.new("L", background.size, 0)
    ImageDraw.Draw(shadow).rounded_rectangle((x, y + 14, x + width, y + height + 14), RADIUS, fill=150)
    shadow = shadow.filter(ImageFilter.GaussianBlur(26))
    base = Image.composite(Image.new("RGB", background.size, (2, 6, 23)), background, shadow)
    mask = Image.new("L", (width * 4, height * 4), 0)
    ImageDraw.Draw(mask).rounded_rectangle((0, 0, width * 4 - 1, height * 4 - 1), RADIUS * 4, fill=255)
    return base, mask.resize((width, height), Image.Resampling.LANCZOS)


def _cursor() -> tuple[Image.Image, tuple[int, int]]:
    scale = 4
    shape = [(0, 0), (0, 21), (5.2, 16.2), (8.6, 24), (12, 22.6), (8.6, 15), (15, 15)]
    pad = 8
    size = (int((15 + pad * 2) * scale), int((25 + pad * 2) * scale))
    layer = Image.new("RGBA", size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(layer)
    points = [((px + pad) * scale, (py + pad) * scale) for px, py in shape]
    shadow = Image.new("RGBA", size, (0, 0, 0, 0))
    ImageDraw.Draw(shadow).polygon([(px + 1.5 * scale, py + 2.5 * scale) for px, py in points], fill=(0, 0, 0, 110))
    layer = Image.alpha_composite(layer, shadow.filter(ImageFilter.GaussianBlur(2.5 * scale)))
    draw = ImageDraw.Draw(layer)
    draw.polygon(points, fill=(255, 255, 255, 255))
    inner = [((px + pad) * scale, (py + pad) * scale) for px, py in
             [(1.6, 3.8), (1.6, 17.2), (5.6, 13.6), (9.2, 21.6), (10.6, 21), (7.2, 13.3), (11.6, 13.3)]]
    draw.polygon(inner, fill=(17, 24, 39, 255))
    final = layer.resize((size[0] // scale * 3 // 2, size[1] // scale * 3 // 2), Image.Resampling.LANCZOS)
    hotspot = (pad * 3 // 2, pad * 3 // 2)
    return final, hotspot


@functools.lru_cache(maxsize=8)
def _white_logo(width: int) -> Image.Image:
    logo = Image.open(LOGO).convert("RGBA")
    alpha = logo.getchannel("A")
    white = Image.new("RGBA", logo.size, (255, 255, 255, 0))
    white.putalpha(alpha)
    height = round(logo.height * width / logo.width)
    return white.resize((width, height), Image.Resampling.LANCZOS)


# -- title-card logo curl ----------------------------------------------------
# The intro card draws the logo's ring in one stroke, as the app's loader does.
# The geometry and the draw-phase easing are read from the app's own markup at
# render time; nothing of it is copied here.

INDEX_HTML = REPO_ROOT / "gbdraw/web/index.html"
CURL_SECONDS = 0.6
_CURL_SUPERSAMPLE = 3
_PATH_TOKENS = re.compile(r"[A-Za-z]|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
Point = tuple[float, float]


def _arc(start: Point, rx: float, ry: float, rotation: float, large: bool, sweep: bool, end: Point) -> list[Point]:
    """Points along an SVG elliptical arc (endpoint parameterization, SVG 1.1 F.6.5)."""

    if not rx or not ry:
        return [end]
    phi = math.radians(rotation)
    cos, sin = math.cos(phi), math.sin(phi)
    dx, dy = (start[0] - end[0]) / 2, (start[1] - end[1]) / 2
    x1, y1 = cos * dx + sin * dy, -sin * dx + cos * dy
    rx, ry = abs(rx), abs(ry)
    scale = x1 * x1 / (rx * rx) + y1 * y1 / (ry * ry)
    if scale > 1:
        rx, ry = rx * math.sqrt(scale), ry * math.sqrt(scale)
    numerator = rx * rx * ry * ry - rx * rx * y1 * y1 - ry * ry * x1 * x1
    coefficient = math.sqrt(max(0.0, numerator / (rx * rx * y1 * y1 + ry * ry * x1 * x1)))
    if large == sweep:
        coefficient = -coefficient
    cx1, cy1 = coefficient * rx * y1 / ry, -coefficient * ry * x1 / rx
    cx = cos * cx1 - sin * cy1 + (start[0] + end[0]) / 2
    cy = sin * cx1 + cos * cy1 + (start[1] + end[1]) / 2
    theta = math.atan2((y1 - cy1) / ry, (x1 - cx1) / rx)
    delta = math.atan2((-y1 - cy1) / ry, (-x1 - cx1) / rx) - theta
    if sweep and delta < 0:
        delta += 2 * math.pi
    elif not sweep and delta > 0:
        delta -= 2 * math.pi
    count = max(2, math.ceil(abs(delta) / math.radians(1.5)))
    points = []
    for step in range(1, count):
        angle = theta + delta * step / count
        x, y = rx * math.cos(angle), ry * math.sin(angle)
        points.append((cos * x - sin * y + cx, sin * x + cos * y + cy))
    return points + [end]


def _path_points(d: str) -> list[list[Point]]:
    """Flatten an SVG path's subpaths (M L H V C A Z, absolute or relative) to points."""

    tokens = _PATH_TOKENS.findall(d)
    subpaths: list[list[Point]] = []
    current = start = (0.0, 0.0)
    command, index = "", 0

    def numbers(count: int) -> list[float]:
        nonlocal index
        values = [float(token) for token in tokens[index:index + count]]
        if len(values) != count or any(token.isalpha() for token in tokens[index:index + count]):
            raise ValueError(f"Path command {command} lacks numbers: {d[:60]}")
        index += count
        return values

    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
        elif not command:
            raise ValueError(f"Path does not start with a command: {d[:60]}")
        relative = command.islower()
        origin = current if relative else (0.0, 0.0)
        kind = command.upper()
        if kind == "Z":
            subpaths[-1].append(start)
            current = start
            command = ""
            continue
        if kind == "M":
            x, y = numbers(2)
            current = start = (origin[0] + x, origin[1] + y)
            subpaths.append([current])
            command = "l" if relative else "L"
            continue
        if not subpaths:
            raise ValueError(f"Path draws before its first move: {d[:60]}")
        if kind == "L":
            x, y = numbers(2)
            points = [(origin[0] + x, origin[1] + y)]
        elif kind == "H":
            (x,) = numbers(1)
            points = [(origin[0] + x, current[1])]
        elif kind == "V":
            (y,) = numbers(1)
            points = [(current[0], origin[1] + y)]
        elif kind == "C":
            values = numbers(6)
            p1, p2, p3 = [(origin[0] + values[i], origin[1] + values[i + 1]) for i in (0, 2, 4)]
            p0 = current
            points = []
            for step in range(1, 49):
                t = step / 48
                a, b, c, e = (1 - t) ** 3, 3 * (1 - t) ** 2 * t, 3 * (1 - t) * t * t, t ** 3
                points.append((a * p0[0] + b * p1[0] + c * p2[0] + e * p3[0],
                               a * p0[1] + b * p1[1] + c * p2[1] + e * p3[1]))
        elif kind == "A":
            rx, ry, rotation, large, sweep, x, y = numbers(7)
            points = _arc(current, rx, ry, rotation, bool(large), bool(sweep), (origin[0] + x, origin[1] + y))
        else:
            raise ValueError(f"Unsupported path command {command}: {d[:60]}")
        subpaths[-1].extend(points)
        current = points[-1]
    if not subpaths:
        raise ValueError("Empty path")
    return subpaths


def _matrix(transform: str | None) -> tuple[float, ...]:
    if not transform:
        return (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    match = re.fullmatch(r"\s*matrix\(([^)]*)\)\s*", transform)
    values = [float(value) for value in re.split(r"[\s,]+", match.group(1).strip())] if match else []
    if len(values) != 6:
        raise ValueError(f"Unsupported transform {transform!r}")
    return tuple(values)


class _TemplateParser(HTMLParser):
    """Collect the loader template's viewBox, pens (per mask), and masked logo paths."""

    def __init__(self) -> None:
        super().__init__()
        self.view_box: str | None = None
        self.pens: dict[str, dict] = {}
        self.shapes: dict[str, list[dict]] = {}
        self._mask: str | None = None
        self._group: str | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if tag == "svg":
            self.view_box = values.get("viewbox")
        elif tag == "mask":
            self._mask = values.get(":id")
        elif tag == "g":
            match = re.search(r"\$\{([\w.]+)\}", values.get(":mask") or "")
            self._group = match.group(1) if match else None
        elif tag == "path" and self._mask:
            self.pens[self._mask] = values
        elif tag == "path" and self._group:
            self.shapes.setdefault(self._group, []).append(values)

    def handle_endtag(self, tag: str) -> None:
        if tag == "mask":
            self._mask = None
        elif tag == "g":
            self._group = None


def _cubic_bezier(x1: float, y1: float, x2: float, y2: float):
    """A CSS ``cubic-bezier()`` timing function."""

    def axis(a: float, b: float, u: float) -> float:
        return 3 * a * u * (1 - u) ** 2 + 3 * b * u * u * (1 - u) + u ** 3

    def timing(x: float) -> float:
        if x <= 0 or x >= 1:
            return min(1.0, max(0.0, x))
        low, high = 0.0, 1.0
        for _ in range(48):
            middle = (low + high) / 2
            low, high = (middle, high) if axis(x1, x2, middle) < x else (low, middle)
        return axis(y1, y2, (low + high) / 2)

    return timing


@dataclass(frozen=True)
class LoaderPen:
    """One mask pen of the loader and the logo part it reveals, in the template's user units."""

    points: tuple[Point, ...]
    distances: tuple[float, ...]   # running length along ``points``
    path_length: float              # the pen's ``pathLength``: the stroke units of its dash interval
    offset: float                   # where the pen starts on the combined stroke, in those units
    shapes: tuple[tuple[Point, ...], ...]

    def drawn(self, stroke: float) -> list[Point]:
        """The pen's points up to ``stroke`` units of the combined stroke."""

        length = min(self.path_length, max(0.0, stroke - self.offset)) * self.distances[-1] / self.path_length
        if length <= 0:
            return []
        cut = bisect.bisect_left(self.distances, length)
        points = list(self.points[:cut])
        if cut < len(self.points):
            (x0, y0), (x1, y1) = self.points[cut - 1], self.points[cut]
            span = self.distances[cut] - self.distances[cut - 1]
            share = (length - self.distances[cut - 1]) / span if span else 1.0
            points.append((x0 + (x1 - x0) * share, y0 + (y1 - y0) * share))
        return points


@dataclass(frozen=True)
class LoaderMark:
    """The app loader's one-stroke logo, parsed from ``#gbdraw-loader-template``."""

    view_box: tuple[float, float, float, float]
    pens: tuple[LoaderPen, ...]
    width: float                    # the pens' stroke width
    total: float                    # the combined stroke's length in stroke units
    easing: tuple[float, float, float, float]

    def drawn(self, seconds: float) -> float:
        """Share of the stroke drawn ``seconds`` into the intro card."""

        return _cubic_bezier(*self.easing)(min(1.0, max(0.0, seconds / CURL_SECONDS)))

    def bounds(self) -> tuple[float, float, float, float]:
        points = [point for pen in self.pens for shape in pen.shapes for point in shape]
        xs, ys = [x for x, _ in points], [y for _, y in points]
        return min(xs), min(ys), max(xs), max(ys)


def loader_mark(html: str | None = None) -> LoaderMark:
    """Parse the loader's pens, logo paths, stroke units, and draw easing from the app's index.html."""

    if html is None:
        return _app_loader_mark()
    template = re.search(r'<script type="text/x-template" id="gbdraw-loader-template">(.*?)</script>', html, re.S)
    easing = re.search(r"@keyframes gbdraw-loader-stroke\s*\{\s*0%\s*\{[^}]*?"
                       r"animation-timing-function:\s*cubic-bezier\(([^)]*)\)", html)
    width = re.search(r"\.gbdraw-loader-pen\s*\{[^}]*?stroke-width:\s*([\d.]+)", html)
    if not (template and easing and width):
        raise ValueError("index.html lacks the loader template, its draw easing, or its pen width")
    parser = _TemplateParser()
    parser.feed(template.group(1))
    view_box = [float(value) for value in (parser.view_box or "").split()]
    if len(view_box) != 4 or not parser.pens or set(parser.pens) != set(parser.shapes):
        raise ValueError("The loader template's viewBox, masks, or masked paths do not parse")
    pens = []
    for key, pen in parser.pens.items():
        name = next((part for part in (pen.get("class") or "").split() if re.fullmatch(r"gbdraw-loader-pen-\w+", part)),
                    None)
        rule = name and re.search(re.escape("." + name) + r"\s*\{\s*--gbdraw-loader-a:\s*clamp\(0,\s*"
                                  r"var\(--gbdraw-loader-from\)(?:\s*-\s*([\d.]+))?,\s*([\d.]+)\)", html)
        if not rule or not pen.get("pathlength") or not pen.get("d"):
            raise ValueError(f"The loader pen {key} has no class rule, pathLength, or path")
        path_length = float(pen["pathlength"])
        if abs(float(rule.group(2)) - path_length) > 1e-6:
            raise ValueError(f"The loader pen {name} clamps to {rule.group(2)}, not its pathLength {path_length}")
        points = tuple(point for subpath in _path_points(pen["d"]) for point in subpath)
        distances = [0.0]
        for (x0, y0), (x1, y1) in zip(points, points[1:]):
            distances.append(distances[-1] + math.hypot(x1 - x0, y1 - y0))
        shapes = []
        for shape in parser.shapes[key]:
            a, b, c, d, e, f = _matrix(shape.get("transform"))
            shapes += [tuple((a * x + c * y + e, b * x + d * y + f) for x, y in subpath)
                       for subpath in _path_points(shape.get("d") or "")]
        pens.append(LoaderPen(points, tuple(distances), path_length, float(rule.group(1) or 0), tuple(shapes)))
    pens.sort(key=lambda pen: pen.offset)
    curve = tuple(float(value) for value in easing.group(1).split(","))
    if len(curve) != 4:
        raise ValueError(f"Unreadable loader easing cubic-bezier({easing.group(1)})")
    mark = LoaderMark(tuple(view_box), tuple(pens), float(width.group(1)),
                      max(pen.offset + pen.path_length for pen in pens), curve)
    left, top, right, bottom = mark.bounds()
    x, y, w, h = mark.view_box
    if left < x or top < y or right > x + w or bottom > y + h:
        raise ValueError("The loader's logo paths leave its viewBox")
    return mark


@functools.lru_cache(maxsize=1)
def _app_loader_mark() -> LoaderMark:
    return loader_mark(INDEX_HTML.read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=1)
def _logo_placement() -> tuple[int, float, float, float]:
    """Where the loader's mark sits on the title logo image.

    Returns the logo's first empty column after the mark (the mark/word split),
    and the scale and offset from loader user units to logo pixels, fitted to
    the mark's bounding box on the image.
    """

    # Half coverage is the edge of a shape, so the antialiased rim does not widen the fit.
    alpha = Image.open(LOGO).getchannel("A").point(lambda value: 255 if value >= 128 else 0)
    columns = list(alpha.resize((alpha.width, 1), Image.Resampling.BOX).getdata())
    first = next(x for x, value in enumerate(columns) if value)
    end = next((x for x in range(first, len(columns)) if not columns[x]), None)
    if end is None:
        raise ValueError(f"{LOGO.name} has no gap between its mark and its word")
    word = next(x for x in range(end, len(columns)) if columns[x])
    box = alpha.crop((first, 0, end, alpha.height)).getbbox()
    left, top, right, bottom = _app_loader_mark().bounds()
    sx, sy = (end - first) / (right - left), (box[3] - box[1]) / (bottom - top)
    if abs(sx / sy - 1) > 0.02:
        raise ValueError(f"The loader's mark does not have the shape of the mark on {LOGO.name}")
    scale = (sx + sy) / 2
    return (end + word) // 2, scale, (first + end) / 2 - scale * (left + right) / 2, \
        (box[1] + box[3]) / 2 - scale * (top + bottom) / 2


def _stroke_outline(points: list[Point], half: float) -> list[Point]:
    """The outline of a butt-capped stroke along a smooth polyline, as one polygon without seams."""

    left, right = [], []
    for index, (x, y) in enumerate(points):
        (x0, y0), (x1, y1) = points[max(0, index - 1)], points[min(len(points) - 1, index + 1)]
        length = math.hypot(x1 - x0, y1 - y0) or 1.0
        nx, ny = -(y1 - y0) / length * half, (x1 - x0) / length * half
        left.append((x + nx, y + ny))
        right.append((x - nx, y - ny))
    return left + right[::-1]


class LogoCurl:
    """The white title logo at one width, split into its word and a mark drawn as the loader draws it."""

    def __init__(self, width: int) -> None:
        self.loader = _app_loader_mark()
        logo = _white_logo(width)
        self.size = logo.size
        split, self._scale, self._dx, self._dy = _logo_placement()
        self._factor = width / Image.open(LOGO).width
        self.split = round(split * self._factor)
        word = logo.getchannel("A")
        word.paste(0, (0, 0, self.split, logo.height))
        self.word = logo.copy()
        self.word.putalpha(word)
        self._marks: dict[float, Image.Image] = {}

    def _point(self, point: Point) -> Point:
        k = self._factor * _CURL_SUPERSAMPLE
        return (point[0] * self._scale + self._dx) * k, (point[1] * self._scale + self._dy) * k

    def mark(self, share: float) -> Image.Image:
        """The mark with ``share`` of its stroke drawn, on the logo's frame."""

        share = round(min(1.0, max(0.0, share)), 4)
        if share in self._marks:
            return self._marks[share]
        big = (self.size[0] * _CURL_SUPERSAMPLE, self.size[1] * _CURL_SUPERSAMPLE)
        alpha = Image.new("L", big, 0)
        stroke = self.loader.width * self._scale * self._factor * _CURL_SUPERSAMPLE
        for pen in self.loader.pens:
            drawn = pen.drawn(share * self.loader.total)
            if len(drawn) < 2:
                continue
            reveal = Image.new("L", big, 0)
            ImageDraw.Draw(reveal).polygon(_stroke_outline([self._point(point) for point in drawn], stroke / 2),
                                           fill=255)
            fill = Image.new("L", big, 0)
            for shape in pen.shapes:
                ImageDraw.Draw(fill).polygon([self._point(point) for point in shape], fill=255)
            alpha = ImageChops.lighter(alpha, ImageChops.multiply(reveal, fill))
        layer = Image.new("RGBA", self.size, (255, 255, 255, 0))
        layer.putalpha(alpha.resize(self.size, Image.Resampling.BOX))
        if share in (0.0, 1.0):
            self._marks[share] = layer
        return layer


@functools.lru_cache(maxsize=4)
def logo_curl(width: int) -> LogoCurl:
    return LogoCurl(width)


class Overlay:
    """Pre-rendered captions, toasts, and badges, blended per frame."""

    def __init__(self, type_: Type, layout: Layout, lang: Language) -> None:
        self.type = type_
        self.layout = layout
        self.lang = lang
        self._cache: dict = {}
        self.brand = self._brand()

    def _brand(self) -> Image.Image:
        text = "gbdraw.app"
        if self.layout.portrait:
            logo = _white_logo(190)
            font = self.type.font(500, 30)
            width = ImageDraw.Draw(Image.new("RGBA", (1, 1))).textlength(text, font=font)
            layer = Image.new("RGBA", (logo.width + 28 + int(width) + 4, max(logo.height, 40)), (0, 0, 0, 0))
            layer.alpha_composite(logo, (0, (layer.height - logo.height) // 2))
            ImageDraw.Draw(layer).text((logo.width + 28, layer.height / 2), text, font=font, anchor="lm",
                                       fill=MUTED + (255,))
            return layer
        logo = _white_logo(150)
        font = self.type.font(500, 22)
        layer = Image.new("RGBA", (logo.width + 4, logo.height + 34), (0, 0, 0, 0))
        layer.alpha_composite(logo, (0, 0))
        draw = ImageDraw.Draw(layer)
        text_width = draw.textlength(text, font=font)
        draw.text(((logo.width - text_width) / 2 + 2, logo.height + 6), text, font=font, fill=MUTED + (255,))
        return layer

    # Intro titles: Latin runs wrap as words; any other character (for example
    # CJK) may break on its own, except closing punctuation, which stays on its line.
    _TOKENS = re.compile(r"[(\[]?[A-Za-z0-9][\w.,:/()%+\-–'’]*\s*|\s+|.")
    # Titles also keep a short full-width parenthetical, such as （烟草）, whole.
    _PHRASES = re.compile(r"（[^（）]{1,8}）|" + _TOKENS.pattern)
    _NO_LINE_START = set("，。、；：？！）」』》〉”’·…")

    def _wrap(self, text: str, font: ImageFont.FreeTypeFont, width: float) -> list[str]:
        draw = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
        lines, line = [], ""
        for token in self._PHRASES.findall(text):
            trial = line + token
            if (line.strip() and token.strip() and token not in self._NO_LINE_START
                    and draw.textlength(trial.rstrip(), font=font) > width):
                lines.append(line.rstrip())
                line = token.lstrip()
            else:
                line = trial
        line = line.rstrip()
        return lines + [line] if line else lines

    def caption_lines(self, step: int | None, title: str, subtitle: str) -> tuple[list[str], list[str]]:
        """The caption's title and subtitle lines; raises ValueError when one needs a forbidden break."""

        layout = self.layout
        width = (layout.card[2] if layout.portrait else 1440) - (76 if step is not None else 0)
        probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
        lines = []
        for text, font in ((title, self.type.font(600, layout.title_size)),
                           (subtitle, self.type.font(400, layout.subtitle_size))):
            try:
                lines.append(wrap_clauses(text, functools.partial(probe.textlength, font=font), width)
                             if text else [])
            except ValueError as error:
                raise ValueError(f"{self.lang.code} {layout.name} caption: {error}") from None
        title_top, subtitle_top, height = layout.caption_rows(len(lines[0]), len(lines[1]))
        title_step, subtitle_step = layout.caption_steps
        bottom = (subtitle_top + subtitle_step * len(lines[1]) if lines[1]
                  else title_top + title_step * len(lines[0]))
        if bottom > height:
            raise ValueError(f"{self.lang.code} {layout.name} caption needs {bottom} px of its {height} px band: "
                             f"{title!r} / {subtitle!r}")
        return lines[0], lines[1]

    def caption(self, step: int | None, title: str, subtitle: str) -> Image.Image:
        key = ("caption", step, title, subtitle)
        if key in self._cache:
            return self._cache[key]
        self.type.check(title + subtitle)
        layout = self.layout
        width = layout.card[2] if layout.portrait else 1440
        indent = 76 if step is not None else 0
        title_font = self.type.font(600, layout.title_size)
        subtitle_font = self.type.font(400, layout.subtitle_size)
        title_lines, subtitle_lines = self.caption_lines(step, title, subtitle)
        title_step, subtitle_step = layout.caption_steps
        y, subtitle_y, height = layout.caption_rows(len(title_lines), len(subtitle_lines))
        layer = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        if step is not None:
            draw.ellipse((0, 24, 56, 80), fill=ACCENT + (255,))
            number = self.type.font(700, 30)
            label = str(step)
            draw.text((28 - draw.textlength(label, font=number) / 2, 32), label, font=number,
                      fill=(255, 255, 255, 255))
        for line in title_lines:
            draw.text((indent, y), line, font=title_font, fill=INK + (255,))
            y += title_step
        y = subtitle_y
        for line in subtitle_lines:
            draw.text((indent, y), line, font=subtitle_font, fill=MUTED + (255,))
            y += subtitle_step
        self._cache[key] = layer
        return layer

    def _scaled(self, layer: Image.Image) -> Image.Image:
        scale = self.layout.ui_scale
        if scale == 1:
            return layer
        return layer.resize((round(layer.width * scale), round(layer.height * scale)), Image.Resampling.LANCZOS)

    def toast(self, text: str, icon: str) -> Image.Image:
        key = ("toast", text, icon)
        if key in self._cache:
            return self._cache[key]
        self.type.check(text)
        probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
        # A long file name takes a smaller type so the toast stays inside the window.
        room = (self.layout.card[2] - 24) / self.layout.ui_scale - 150
        size = 26
        while size > 16 and probe.textlength(text, font=self.type.font(500, size)) > room:
            size -= 1
        font = self.type.font(500, size)
        width = int(probe.textlength(text, font=font)) + 110
        scale = 3
        layer = Image.new("RGBA", ((width + 40) * scale, 108 * scale), (0, 0, 0, 0))
        shadow = Image.new("L", layer.size, 0)
        ImageDraw.Draw(shadow).rounded_rectangle((20 * scale, 26 * scale, (width + 20) * scale, 94 * scale),
                                                  34 * scale, fill=90)
        layer.putalpha(shadow.filter(ImageFilter.GaussianBlur(10 * scale)))
        layer = Image.merge("RGBA", (*Image.new("RGB", layer.size, (2, 6, 23)).split(), layer.getchannel("A")))
        draw = ImageDraw.Draw(layer)
        draw.rounded_rectangle((20 * scale, 18 * scale, (width + 20) * scale, 86 * scale), 34 * scale,
                               fill=(255, 255, 255, 250), outline=(203, 213, 225, 255), width=scale * 2)
        cx, cy = 60 * scale, 52 * scale
        draw.ellipse((cx - 17 * scale, cy - 17 * scale, cx + 17 * scale, cy + 17 * scale), fill=(22, 163, 74, 255))
        if icon == "download":
            draw.line([(cx, cy - 9 * scale), (cx, cy + 6 * scale)], fill="white", width=4 * scale)
            draw.line([(cx - 7 * scale, cy - 1 * scale), (cx, cy + 7 * scale), (cx + 7 * scale, cy - 1 * scale)],
                      fill="white", width=4 * scale, joint="curve")
        else:
            draw.line([(cx - 8 * scale, cy), (cx - 2 * scale, cy + 6 * scale), (cx + 9 * scale, cy - 6 * scale)],
                      fill="white", width=4 * scale, joint="curve")
        layer = layer.resize((width + 40, 108), Image.Resampling.LANCZOS)
        ImageDraw.Draw(layer).text((90, 33 + (26 - size) / 2), text, font=font, fill=(15, 23, 42, 255))
        self._cache[key] = self._scaled(layer)
        return self._cache[key]

    def badge(self, factor: float) -> Image.Image:
        text = f"{factor:.0f}\u00d7" if factor >= 2 else self.lang("fast")
        key = ("badge", text)
        if key in self._cache:
            return self._cache[key]
        font = self.type.font(600, 24)
        probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
        width = int(probe.textlength(text, font=font)) + 78
        scale = 3
        layer = Image.new("RGBA", (width * scale, 48 * scale), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        draw.rounded_rectangle((0, 0, width * scale - 1, 48 * scale - 1), 24 * scale, fill=(15, 23, 42, 215))
        for offset in (0, 13):
            x = (20 + offset) * scale
            draw.polygon([(x, 14 * scale), (x, 34 * scale), (x + 13 * scale, 24 * scale)], fill=(255, 255, 255, 255))
        layer = layer.resize((width, 48), Image.Resampling.LANCZOS)
        ImageDraw.Draw(layer).text((56, 9), text, font=font, fill=(255, 255, 255, 255))
        self._cache[key] = self._scaled(layer)
        return self._cache[key]


_ITALIC = re.compile(r"<i>(.*?)</i>")


def _runs(text: str) -> list[tuple[str, bool]]:
    """Split ``<i>`` markup (Gallery titles italicize taxa) into styled runs."""

    runs, last = [], 0
    for match in _ITALIC.finditer(text):
        runs += [(text[last:match.start()], False), (match.group(1), True)]
        last = match.end()
    runs.append((text[last:], False))
    return [run for run in runs if run[0]]


# Chinese has no spaces, so a break between two Han characters in an intro
# title is kept only after punctuation or the particle 的, where it cannot
# split a word. Captions are stricter (``wrap_clauses``).
_SOFT_BREAK = set("，。、；：？！）」』》〉”’的")

# A caption line ends only at a clause or sentence boundary (Owner rule of
# 2026-10-11, VR-D09): after ". ; : , ? !" and a space, after "· › —", or after
# their CJK forms. A word or a clause is never split, and 的 is no boundary.
_CLAUSE_BREAK = re.compile(r"(?<=[.;:,?!])\s+|(?<=[·›—])\s*|(?<=[，。；：、？！])(?![」』》〉”’）])\s*")


def wrap_clauses(text: str, measure, width: float) -> list[str]:
    """Fill lines of at most ``width`` with whole clauses of ``text``.

    ``measure`` returns a string's width. A clause wider than ``width`` raises
    ValueError naming the caption, so its author shortens it.
    """

    cuts = sorted({match.end() for match in _CLAUSE_BREAK.finditer(text)} - {0, len(text)})
    edges = [0, *cuts, len(text)]
    lines, line = [], ""
    for clause in (text[a:b] for a, b in zip(edges, edges[1:])):
        if line and measure((line + clause).rstrip()) > width:
            lines.append(line.rstrip())
            line = clause
        else:
            line += clause
    lines.append(line.rstrip())
    for line in lines:
        if measure(line) > width:
            raise ValueError(f"{text!r} needs a line break inside {line!r} ({measure(line):.0f} px, "
                             f"room for {width:.0f} px); shorten it, or end a clause earlier with "
                             ". ; : , · › — or their CJK forms")
    return lines


def _clean_breaks(plain: str, lines: list[str]) -> bool:
    cursor = 0
    for line in lines[:-1]:
        cursor = plain.index(line, cursor) + len(line)
        before, after = plain[cursor - 1], plain[cursor]
        if not (after.isspace() or after == "（" or before in _SOFT_BREAK or before.isascii() or after.isascii()):
            return False
    return True


def _fade(layer: Image.Image, alpha: float) -> Image.Image:
    if alpha >= 0.999:
        return layer
    faded = layer.copy()
    faded.putalpha(layer.getchannel("A").point(lambda v: round(v * max(0.0, alpha))))
    return faded




def _menu(draw: ImageDraw.ImageDraw, type_: Type, event: dict, to_screen, scale: float, alpha: float) -> None:
    """Draw the native option list that headless Chromium leaves unpainted."""

    left, top = to_screen(event["x"], event["y"] + event["height"] + 2)
    first, rows = event["window"]
    labels = event["options"][first:first + rows]
    font = type_.font(400, max(9, round(13 * scale)))
    row = 22 * scale
    width = max(event["width"] * scale,
                max(draw.textlength(label, font=font) for label in labels) + 28 * scale)
    height = row * len(labels) + 8 * scale
    a = round(255 * alpha)
    draw.rounded_rectangle((left + 2, top + 5, left + width + 2, top + height + 5), 7 * scale,
                           fill=(15, 23, 42, round(40 * alpha)))
    draw.rounded_rectangle((left, top, left + width, top + height), 7 * scale,
                           fill=(255, 255, 255, a), outline=(203, 213, 225, a), width=max(1, round(scale)))
    for index, label in enumerate(labels):
        y = top + 4 * scale + index * row
        chosen = first + index == event["selected"]
        if chosen:
            draw.rounded_rectangle((left + 4 * scale, y, left + width - 4 * scale, y + row - 2 * scale),
                                   4 * scale, fill=(37, 99, 235, a))
        draw.text((left + 12 * scale, y + row / 2 - 1 * scale), label, font=type_.label_font(label, font.size),
                  anchor="lm",
                  fill=(255, 255, 255, a) if chosen else (30, 41, 59, a))


def _swatch(draw: ImageDraw.ImageDraw, type_: Type, event: dict, to_screen, scale: float, alpha: float) -> None:
    left, top = to_screen(event["x"] - 70, event["y"] + event["height"] + 4)
    width, height = 176 * scale, 58 * scale
    a = round(255 * alpha)
    color = tuple(int(event["value"][i:i + 2], 16) for i in (1, 3, 5))
    draw.rounded_rectangle((left + 2, top + 5, left + width + 2, top + height + 5), 8 * scale,
                           fill=(15, 23, 42, round(40 * alpha)))
    draw.rounded_rectangle((left, top, left + width, top + height), 8 * scale,
                           fill=(255, 255, 255, a), outline=(203, 213, 225, a), width=max(1, round(scale)))
    draw.rounded_rectangle((left + 9 * scale, top + 9 * scale, left + 49 * scale, top + 49 * scale), 6 * scale,
                           fill=color + (a,))
    draw.text((left + 60 * scale, top + height / 2), event["value"].upper(),
              font=type_.font(500, max(9, round(15 * scale))), anchor="lm", fill=(30, 41, 59, a))


# -- recording --------------------------------------------------------------

class Recording:
    def __init__(self, raw: Path) -> None:
        self.raw = raw
        self.meta = json.loads((raw / "walkthrough.json").read_text(encoding="utf-8"))
        self.story = STORIES[self.meta["journey"]]
        self.frames = json.loads((raw / "frames.json").read_text(encoding="utf-8"))
        self.events = json.loads((raw / "events.json").read_text(encoding="utf-8"))
        self.viewport = tuple(self.meta["viewport"])
        self.scale = self.meta["scale"]
        begin = self.frames[0][0]
        end = next(e["t"] for e in self.events if e["kind"] == "end")
        self.timemap = TimeMap.build(self.events, begin, end)
        self.pointer = Pointer(self.events, self.timemap)
        self.frame_times = [row[0] for row in self.frames]
        out = self.timemap.out
        self.captions = [(out(e["t"]), (e["step"], e["title"], e["subtitle"]))
                         for e in self.events if e["kind"] == "caption"]
        self.caption_times = [c[0] for c in self.captions]
        self.toasts = [(out(e["t"]), (e["text"], e["icon"])) for e in self.events if e["kind"] == "toast"]
        self.popovers = [(out(e["t"]), out(e["t"] + e["duration"]), e)
                         for e in self.events if e["kind"] in ("menu", "swatch")]
        self.marks = {e["name"]: out(e["t"]) for e in self.events if e["kind"] == "mark"}
        self.steps = sorted({step for _, (step, _, _) in self.captions if step is not None})
        self._decoded: tuple[int, Image.Image] | None = None

    @property
    def duration(self) -> float:
        return self.timemap.duration

    def source_index(self, moment: float) -> int:
        return max(0, bisect.bisect_right(self.frame_times, self.timemap.source(moment)) - 1)

    def source_frame(self, moment: float) -> Image.Image:
        index = self.source_index(moment)
        if self._decoded is None or self._decoded[0] != index:
            with Image.open(self.raw / "frames" / self.frames[index][1]) as image:
                self._decoded = (index, image.convert("RGB"))
        return self._decoded[1]

    def caption_at(self, moment: float, since: float) -> tuple[tuple, float, tuple | None]:
        """Caption at ``moment``, its fade-in progress, and the caption it replaces.

        ``since`` is where the current cut began; a caption that was already
        showing then simply fades in with the cut.
        """

        index = bisect.bisect_right(self.caption_times, moment) - 1
        if index < 0:
            return self.captions[0][1], 1.0, None
        start, current = self.captions[index]
        previous = self.captions[index - 1][1] if index and start > since else None
        return current, (moment - max(start, since)) / 0.35, previous


# -- edit plans -------------------------------------------------------------

@dataclass
class Piece:
    kind: str                 # "intro", "montage", "recording", "finale", or "outro"
    frames: int
    moments: tuple = ()       # recording pieces: walkthrough-timeline position of each frame
    speed: float = 1.0
    fade: int = 0             # frames of cross-dissolve from the previous piece
    caption: tuple | None = None   # recording pieces: a fixed caption instead of the recording's

    @property
    def start(self) -> float:
        return self.moments[0] if self.moments else 0.0

    def moment(self, local: int) -> float:
        if local < len(self.moments):
            return self.moments[local]
        return self.moments[-1] + (local - len(self.moments) + 1) * self.speed / FPS

    def rate(self, local: int) -> float:
        return (self.moment(local + 1) - self.moment(local)) * FPS


def _moments(timemap: TimeMap, start: float, stop: float, speed: float, wait_speed: float) -> tuple:
    moments, moment = [], start
    while moment < stop:
        moments.append(moment)
        moment += speed * (wait_speed if timemap.waiting(moment) else 1.0) / FPS
    return tuple(moments)


@dataclass
class Edit:
    name: str
    pieces: list[Piece]
    numbered: bool            # show step numbers and the progress bar
    badge_from: float         # show the fast-forward badge at or above this factor
    finale_caption: tuple


# Highlights cut the walkthrough by named marks, so a new recording with other
# timings still yields the same story. Offsets and speeds are on the
# walkthrough timeline.
HIGHLIGHTS = (
    Clip("upload", -0.3, "loaded", 1.2),
    Clip("generate", 0.0, "generated", 2.0),
    Clip("layout", -0.2, "layout-done", 2.0),
    Clip("labels", -0.2, "crowded", 2.2),
    Clip("priority-type", -0.2, "genes", 2.6),
    Clip("rule-4", -0.1, "legend", 2.4),
    Clip("bracket", -0.2, "dloop", 2.6),
    Clip("export", -0.1, "downloaded", 1.4),
)
HIGHLIGHT_SPEED = 1.4
HIGHLIGHT_WAIT_SPEED = 2.5  # generation waits, on top of the walkthrough's own compression
MONTAGE_SECONDS = 3.6
HIGHLIGHT_OUTRO_SECONDS = 4.2
CUT_FADE = 7
CHAPTER_FADE = round(FADE_SECONDS * FPS)


def _clip_pieces(recording: Recording, clips: tuple[Clip, ...], speed: float, wait_speed: float,
                 first_fade: int) -> list[Piece]:
    """One recording piece per clip, cut at the recording's marks; later clips dissolve in quickly."""

    pieces = []
    for index, clip in enumerate(clips):
        missing = {clip.first, clip.last} - set(recording.marks)
        if missing:
            raise ValueError(f"Recording lacks marks {sorted(missing)}")
        start = recording.marks[clip.first] + clip.before
        stop = recording.marks[clip.last] + clip.after
        if stop <= start:
            raise ValueError(f"Clip {clip.first}->{clip.last} is empty")
        moments = _moments(recording.timemap, start, stop, speed, wait_speed)
        pieces.append(Piece("recording", len(moments), moments=moments, speed=speed,
                            fade=first_fade if index == 0 else CUT_FADE,
                            caption=(None, clip.label, "") if clip.label else None))
    return pieces


def full_edit(recording: Recording, finale_frames: int) -> Edit:
    # The vector zoom continues the last step, which exports the SVG. Journeys
    # that end in the exported file itself, such as an Interactive SVG opened
    # in the browser, have no zoom.
    finale_caption = (recording.steps[-1], *FINALE_TEXT)
    opening = _clip_pieces(recording, recording.story.opening, 1.0, 1.0, 0)
    pieces = [
        *opening,
        # After an opening the title card cuts in quickly, so its logo curl is not drawn under the dissolve.
        Piece("intro", round((OPENING_INTRO_SECONDS if opening else INTRO_SECONDS) * FPS),
              fade=CUT_FADE if opening else 0),
        Piece("recording", math.ceil(recording.duration * FPS),
              moments=tuple(n / FPS for n in range(math.ceil(recording.duration * FPS))), fade=CHAPTER_FADE),
    ]
    if finale_frames:
        pieces.append(Piece("finale", finale_frames, fade=CHAPTER_FADE))
    pieces.append(Piece("outro", round(OUTRO_SECONDS * FPS), fade=CHAPTER_FADE))
    return Edit("walkthrough", pieces, numbered=True, badge_from=1.3, finale_caption=finale_caption)


def highlights_edit(recording: Recording, finale_frames: int) -> Edit:
    pieces = [Piece("montage", round(MONTAGE_SECONDS * FPS)),
              *_clip_pieces(recording, HIGHLIGHTS, HIGHLIGHT_SPEED, HIGHLIGHT_WAIT_SPEED, CHAPTER_FADE)]
    pieces.append(Piece("finale", round(finale_frames / 1.2), speed=1.2, fade=CHAPTER_FADE))
    pieces.append(Piece("outro", round(HIGHLIGHT_OUTRO_SECONDS * FPS), fade=CHAPTER_FADE))
    return Edit("highlights", pieces, numbered=False, badge_from=2.0,
                finale_caption=(None, *FINALE_TEXT))


EDITS = {"walkthrough": full_edit, "highlights": highlights_edit}


# -- composition ------------------------------------------------------------

class Composer:
    def __init__(self, recording: Recording, type_: Type, layout: Layout, lang: Language,
                 edits: set[str]) -> None:
        raw = recording.raw
        self.recording = recording
        self.layout = layout
        self.lang = lang
        self.type = type_
        self.overlay = Overlay(type_, layout, lang)
        self.camera = Camera(recording.events, recording.timemap, recording.viewport, layout)
        self._thumbnail: tuple[int, Image.Image] | None = None
        self.background = _background(layout.size)
        self.base, self.mask = _card_layers(self.background, layout.card)
        self.cursor, self.hotspot = _cursor()
        self.finale = sorted((raw / layout.finale).glob("*.png"))
        if not self.finale and layout.finale in recording.meta["finale"]:
            raise FileNotFoundError(f"Missing {layout.finale} frames in {raw}")
        self.figure = Image.open(raw / "final-figure.png").convert("RGB")
        self.gallery = [(Image.open(raw / row["path"]).convert("RGB"), row["label"])
                        for row in recording.meta.get("gallery", [])]
        self.edits = {name: EDITS[name](recording, len(self.finale)) for name in edits}
        self._static: dict = {}
        # Lay out every caption first: one that needs a forbidden line break
        # fails the render before anything is encoded.
        for edit in self.edits.values():
            shown = [caption for _, caption in recording.captions]
            shown += [edit.finale_caption] if any(piece.kind == "finale" for piece in edit.pieces) else []
            shown += [piece.caption for piece in edit.pieces if piece.caption]
            for step, title, subtitle in shown + ([MONTAGE_CAPTION] if edit.name == "highlights" else []):
                self.overlay.caption(step if edit.numbered else None, lang(title), lang(subtitle))

    # recording ------------------------------------------------------------
    def _screen(self, moment: float, view: tuple[float, float, float]) -> Image.Image:
        rec = self.recording
        frame = rec.source_frame(moment)
        cx, cy, vh = view
        vw = vh * self.layout.aspect
        left, top = cx - vw / 2, cy - vh / 2
        s = rec.scale
        width, height = self.layout.card[2], self.layout.card[3]
        screen = frame.resize((width, height), Image.Resampling.BICUBIC,
                              box=(left * s, top * s, (left + vw) * s, (top + vh) * s))
        scale = width / vw

        def to_screen(x: float, y: float) -> tuple[float, float]:
            return (x - left) * scale, (y - top) * scale

        self._inset(screen, moment, view)
        draw = ImageDraw.Draw(screen, "RGBA")
        for start, stop, event in rec.popovers:
            if start <= moment < stop + 0.15:
                alpha = min(_ease((moment - start) / 0.12), _ease((stop + 0.15 - moment) / 0.15))
                (_menu if event["kind"] == "menu" else _swatch)(draw, self.type, event, to_screen, scale, alpha)
        for rx, ry, p in rec.pointer.ripples(moment):
            ox, oy = to_screen(rx, ry)
            radius = 10 + 26 * _ease(p)
            draw.ellipse((ox - radius, oy - radius, ox + radius, oy + radius),
                         outline=ACCENT + (round(150 * (1 - p)),), width=4)
        cursor = self.cursor
        if rec.pointer.pressed(moment):
            cursor = cursor.resize((round(cursor.width * 0.88), round(cursor.height * 0.88)), Image.Resampling.LANCZOS)
        sx, sy = to_screen(*rec.pointer.at(moment))
        hx = self.hotspot[0] * cursor.width / self.cursor.width
        hy = self.hotspot[1] * cursor.height / self.cursor.height
        screen.paste(cursor, (round(sx - hx), round(sy - hy)), cursor)
        return screen

    def _inset_chrome(self, size: tuple[int, int]) -> tuple[Image.Image, Image.Image, Image.Image, int]:
        """The inset's shadow, rounded mask, and border, made once per size."""

        key = ("inset", size)
        if key not in self._static:
            width, height = size
            pad, radius, ss = 16, 9, 4
            outer = (width + 2 * pad, height + 2 * pad)
            shadow = Image.new("L", outer, 0)
            ImageDraw.Draw(shadow).rounded_rectangle((pad, pad + 5, pad + width, pad + height + 5), radius, fill=150)
            shadow_layer = Image.new("RGBA", outer, (2, 6, 23, 0))
            shadow_layer.putalpha(shadow.filter(ImageFilter.GaussianBlur(8)))
            mask = Image.new("L", (width * ss, height * ss), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, width * ss - 1, height * ss - 1), radius * ss, fill=255)
            border = Image.new("RGBA", (outer[0] * ss, outer[1] * ss), (0, 0, 0, 0))
            ImageDraw.Draw(border).rounded_rectangle(
                (pad * ss, pad * ss, (pad + width) * ss - 1, (pad + height) * ss - 1), radius * ss,
                outline=(226, 232, 240, 255), width=2 * ss)
            self._static[key] = (shadow_layer, mask.resize(size, Image.Resampling.LANCZOS),
                                 border.resize(outer, Image.Resampling.LANCZOS), pad)
        return self._static[key]

    def _inset(self, screen: Image.Image, moment: float, view: tuple[float, float, float]) -> None:
        """Draw the whole stage, small, with the part the window shows marked."""

        corner, alpha = self.camera.inset_at(moment)
        if corner is None or alpha <= 0:
            return
        rec = self.recording
        x, y, width, height = self.layout.inset_box(rec.viewport, corner)
        index = rec.source_index(moment)
        if self._thumbnail is None or self._thumbnail[0] != index:
            stage = rec.source_frame(moment).reduce(2)
            self._thumbnail = (index, stage.resize((width, height), Image.Resampling.LANCZOS))
        thumbnail = self._thumbnail[1].convert("RGBA")
        cx, cy, vh = view
        vw, scale = vh * self.layout.aspect, width / rec.viewport[0]
        visible = ((cx - vw / 2) * scale, (cy - vh / 2) * scale, (cx + vw / 2) * scale - 1, (cy + vh / 2) * scale - 1)
        shade = Image.new("RGBA", thumbnail.size, (15, 23, 42, 85))
        ImageDraw.Draw(shade).rectangle(visible, fill=(0, 0, 0, 0))
        thumbnail.alpha_composite(shade)
        ImageDraw.Draw(thumbnail).rectangle(visible, outline=ACCENT + (255,), width=3)
        shadow, mask, border, pad = self._inset_chrome((width, height))
        thumbnail.putalpha(mask)
        area = (x - pad, y - pad, x + width + pad, y + height + pad)
        under = screen.crop(area).convert("RGBA")
        inset = Image.alpha_composite(under, shadow)
        inset.alpha_composite(thumbnail, (pad, pad))
        inset.alpha_composite(border)
        screen.paste(Image.blend(under, inset, alpha).convert("RGB") if alpha < 1 else inset.convert("RGB"), area[:2])

    def _progress(self, canvas: Image.Image, step: int | None) -> None:
        steps = self.recording.steps
        if step is None or not steps:
            return
        draw = ImageDraw.Draw(canvas)
        card = self.layout.card
        # Tall feeds cover the bottom of the frame, so the bar sits under the window there.
        gap = 8
        y = card[1] + card[3] + PROGRESS_GAP if self.layout.portrait else self.layout.size[1] - 14
        width = (card[2] - gap * (len(steps) - 1)) / len(steps)
        for index, value in enumerate(steps):
            x = card[0] + index * (width + gap)
            color = ACCENT if value == step else (INK if value < step else (51, 65, 85))
            draw.rounded_rectangle((x, y, x + width, y + 4), 2, fill=color)

    def _compose_card(self, screen: Image.Image, caption: tuple, caption_in: float, previous: tuple | None,
                      extras: list, edit: Edit) -> Image.Image:
        card = self.layout.card
        canvas = self.base.copy()
        canvas.paste(screen, card[:2], self.mask)
        canvas = canvas.convert("RGBA")
        for layer, position in extras:
            canvas.alpha_composite(layer, position)
        if not edit.numbered:
            caption = (None, *caption[1:])
            previous = previous and (None, *previous[1:])
        if previous is not None and caption_in < 1:
            self._caption(canvas, previous, 1 - _ease(caption_in * 2), 0)
        p = _ease(caption_in)
        if p > 0:
            self._caption(canvas, caption, p, round((1 - p) * 14))
        self._brand(canvas, below_progress=edit.numbered)
        if edit.numbered:
            self._progress(canvas, caption[0])
        return canvas.convert("RGB")

    def _caption(self, canvas: Image.Image, caption: tuple, alpha: float, rise: int) -> None:
        card = self.layout.card
        step, title, subtitle = caption
        layer = self.overlay.caption(step, self.lang(title), self.lang(subtitle))
        if self.layout.portrait:
            # Bottom-aligned just above the window, so one- and two-line titles sit alike.
            position = (card[0], card[1] - 26 - layer.height + rise)
        else:
            position = (card[0], card[1] + card[3] + 4 + rise)
        canvas.alpha_composite(_fade(layer, alpha), position)

    def _brand(self, canvas: Image.Image, *, below_progress: bool = False) -> None:
        card, brand = self.layout.card, self.overlay.brand
        if self.layout.portrait:
            y = card[1] + card[3] + (PROGRESS_GAP + 38 if below_progress else 34)
            canvas.alpha_composite(brand, ((self.layout.size[0] - brand.width) // 2, y))
        else:
            canvas.alpha_composite(brand, (card[0] + card[2] - brand.width, card[1] + card[3] + 34))

    def _recording_frame(self, piece: Piece, local: int, edit: Edit) -> Image.Image:
        rec = self.recording
        moment = piece.moment(local)
        extras = []
        inset = self.camera.inset_at(moment)
        for start, (text, icon) in rec.toasts:
            age = moment - start
            if 0 <= age < 2.8:
                alpha = min(_ease(age / 0.3), _ease((2.8 - age) / 0.4))
                layer = self.overlay.toast(self.lang(text), icon)
                x, y = self.layout.toast_origin(rec.viewport, layer.size, *inset)
                extras.append((_fade(layer, alpha), (x, y + round((1 - alpha) * 16))))
        factor = rec.timemap.factor(moment) * piece.rate(local)
        # An opening clip is a montage of results: no badge, even where it dissolves into a skimmed moment.
        if factor >= edit.badge_from and piece.caption is None:
            badge = self.overlay.badge(factor)
            extras.append((badge, self.layout.badge_origin(badge.size)))
        if piece.caption:
            caption, caption_in, previous = piece.caption, local / FPS / 0.35, None
        else:
            caption, caption_in, previous = rec.caption_at(moment, piece.start)
        return self._compose_card(self._screen(moment, self.camera.at(moment)), caption, caption_in, previous,
                                  extras, edit)

    def _finale_frame(self, piece: Piece, local: int, edit: Edit) -> Image.Image:
        index = min(len(self.finale) - 1, round(local * piece.speed))
        with Image.open(self.finale[index]) as image:
            screen = image.convert("RGB").resize(self.layout.card[2:], Image.Resampling.LANCZOS)
        return self._compose_card(screen, edit.finale_caption, local / FPS / 0.35, None, [], edit)

    # cards ----------------------------------------------------------------
    def _text(self, canvas: Image.Image, text: str, weight: int, size: int, color: tuple, alpha: float,
              position: tuple[float, float], *, center: bool = False, rise: int = 12) -> None:
        if alpha <= 0:
            return
        runs = _runs(text)
        self.type.check("".join(run for run, _ in runs))
        font = self.type.font(weight, size)
        key = ("text", text, weight, size, color)
        if key not in self._static:
            probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
            height = size + 24
            lean = round(SLANT * height)
            layer = Image.new("RGBA", (int(sum(probe.textlength(run, font=font) for run, _ in runs)) + 4 + lean,
                                       height), (0, 0, 0, 0))
            x = 0.0
            for run, italic in runs:
                piece = Image.new("RGBA", (int(probe.textlength(run, font=font)) + 4 + lean, height), (0, 0, 0, 0))
                ImageDraw.Draw(piece).text((2, 0), run, font=font, fill=color + (255,))
                if italic:
                    # Lean the run to the right about its baseline.
                    baseline = font.getmetrics()[0]
                    piece = piece.transform(piece.size, Image.Transform.AFFINE,
                                            (1, SLANT, -SLANT * baseline, 0, 1, 0), Image.Resampling.BICUBIC)
                layer.alpha_composite(piece, (round(x), 0))
                x += probe.textlength(run, font=font)
            self._static[key] = layer
        layer = self._static[key]
        x, y = position
        if center:
            x -= layer.width / 2
        canvas.alpha_composite(_fade(layer, alpha), (round(x), round(y + (1 - alpha) * rise)))

    def _intro_frame(self, local: int) -> Image.Image:
        t = local / FPS
        canvas = self.background.copy().convert("RGBA")
        title, note = self.recording.story.intro
        if self.layout.portrait:
            return self._portrait_intro(canvas, t, self.lang(title), self.lang(note))
        middle = self.layout.size[0] / 2
        self._intro_logo(canvas, 820, middle, 300, t)
        # A long title wraps onto more lines instead of running off the frame.
        y = 580
        for line in self._wrapped(self.lang(title), 600, 46, self.layout.size[0] - 320):
            self._text(canvas, line, 600, 46, INK, _ease((t - 0.6) / 0.6), (middle, y), center=True)
            y += 60
        self._text(canvas, self.lang(note), 400, 28, MUTED, _ease((t - 0.9) / 0.6), (middle, y - 60 + 76), center=True)
        return canvas.convert("RGB")

    def _intro_logo(self, canvas: Image.Image, width: int, middle: float, top: int, t: float) -> None:
        """The title logo: its ring drawn in one stroke as the app's loader draws it, then its word."""

        curl = logo_curl(width)
        origin = (round(middle - width / 2), top)
        canvas.alpha_composite(curl.mark(curl.loader.drawn(t)), origin)
        canvas.alpha_composite(_fade(curl.word, _ease((t - 0.45) / 0.45)), origin)

    def _wrapped(self, text: str, weight: int, size: int, width: float) -> list[str]:
        """Wrap ``text`` to ``width`` in even lines, keeping its ``<i>`` runs on each line they reach."""

        runs = _runs(text)
        plain = "".join(run for run, _ in runs)
        slanted = [italic for run, italic in runs for _ in run]
        font = self.type.font(weight, size)
        wrap = functools.partial(self.overlay._wrap, font=font)
        sentences = [part.strip() for part in re.split(r"(?<=[。！？])|(?<=[.!?])\s+", plain) if part.strip()]
        if len(sentences) > 1 and all(len(wrap(part, width=width)) == 1 for part in sentences):
            wrapped = sentences
        else:
            # Narrow the measure while the line count holds, so no line is left
            # nearly empty, and keep the narrowest wrap that splits no word.
            candidates = [wrap(plain, width=width)]
            while len(candidates[-1]) > 1:
                width -= size / 2
                trial = wrap(plain, width=width)
                if len(trial) != len(candidates[0]):
                    break
                candidates.append(trial)
            wrapped = ([lines for lines in candidates if _clean_breaks(plain, lines)] or candidates)[-1]
        lines, cursor = [], 0
        for line in wrapped:
            start = plain.index(line, cursor)
            cursor = start + len(line)
            markup, italic = "", False
            for char, flag in zip(line, slanted[start:cursor]):
                if flag != italic:
                    markup += "<i>" if flag else "</i>"
                    italic = flag
                markup += char
            lines.append(markup + ("</i>" if italic else ""))
        return lines

    def _portrait_intro(self, canvas: Image.Image, t: float, title: str, note: str) -> Image.Image:
        """Logo, wrapped title, and note, centered on the window area of a tall frame."""

        card = self.layout.card
        middle = self.layout.size[0] / 2
        logo = _white_logo(760)
        title_size, note_size = self.layout.title_size, self.layout.subtitle_size
        title_lines = self._wrapped(title, 600, title_size, card[2])
        note_lines = self._wrapped(note, 400, note_size, card[2])
        title_step, note_step = round(title_size * 1.3), round(note_size * 1.5)
        height = logo.height + 110 + title_step * len(title_lines) + 30 + note_step * len(note_lines)
        y = card[1] + (card[3] - height) / 2
        self._intro_logo(canvas, logo.width, middle, round(y), t)
        y += logo.height + 110
        for line in title_lines:
            self._text(canvas, line, 600, title_size, INK, _ease((t - 0.6) / 0.6), (middle, y), center=True)
            y += title_step
        y += 30
        for line in note_lines:
            self._text(canvas, line, 400, note_size, MUTED, _ease((t - 0.9) / 0.6), (middle, y), center=True)
            y += note_step
        return canvas.convert("RGB")

    def _gallery_card(self, index: int, size: tuple[int, int]) -> Image.Image:
        key = ("gallery", index, size)
        if key not in self._static:
            figure, label = self.gallery[index]
            label = self.lang(label)
            width, height = size
            card = Image.new("RGBA", size, (255, 255, 255, 255))
            fitted = figure.copy()
            label_size = 30 if self.layout.portrait else 24
            fitted.thumbnail((width - 36, height - 34 - label_size * 2), Image.Resampling.LANCZOS)
            card.paste(fitted, ((width - fitted.width) // 2,
                                18 + (height - 34 - label_size * 2 - fitted.height) // 2))
            draw = ImageDraw.Draw(card)
            self.type.check(label)
            draw.text((22, height - 22 - label_size), label, font=self.type.font(600, label_size),
                      fill=(15, 23, 42, 255))
            mask = Image.new("L", (width * 3, height * 3), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, width * 3 - 1, height * 3 - 1), 14 * 3, fill=255)
            card.putalpha(mask.resize(size, Image.Resampling.LANCZOS))
            self._static[key] = card
        return self._static[key]

    def _montage_frame(self, local: int) -> Image.Image:
        t = local / FPS
        canvas = self.background.copy().convert("RGBA")
        card = self.layout.card
        columns, rows = (2, 3) if self.layout.portrait else (3, 2)
        gap = 24 if self.layout.portrait else 28
        width = (card[2] - gap * (columns - 1)) // columns
        height = (card[3] - gap * (rows - 1)) // rows
        for index in range(min(len(self.gallery), columns * rows)):
            p = _ease((t + 0.25 - index * 0.14) / 0.45)
            if p <= 0:
                continue
            tile = self._gallery_card(index, (width, height))
            grow = 0.94 + 0.06 * p
            scaled = tile.resize((round(width * grow), round(height * grow)), Image.Resampling.BILINEAR)
            x = card[0] + (index % columns) * (width + gap) + (width - scaled.width) // 2
            y = card[1] + (index // columns) * (height + gap) + (height - scaled.height) // 2
            canvas.alpha_composite(_fade(scaled, p), (x, y))
        self._caption(canvas, MONTAGE_CAPTION, _ease((t - 0.8) / 0.5), 0)
        self._brand(canvas)
        return canvas.convert("RGB")

    def _outro_card(self, box: tuple[int, int]) -> Image.Image:
        key = ("outro-card", box)
        if key not in self._static:
            figure = self.figure.copy()
            figure.thumbnail((box[0] - 56, box[1] - 56), Image.Resampling.LANCZOS)
            card = Image.new("RGBA", (figure.width + 56, figure.height + 56), (255, 255, 255, 255))
            card.paste(figure, (28, 28))
            mask = Image.new("L", (card.width * 3, card.height * 3), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, mask.width - 1, mask.height - 1), RADIUS * 3, fill=255)
            card.putalpha(mask.resize(card.size, Image.Resampling.LANCZOS))
            self._static[key] = card
        return self._static[key]

    OUTRO_ROWS = (("gbdraw.app", 700, 60, INK, 0.6),
                  ("Free and open source (MIT)", 400, 32, MUTED, 0.9),
                  ("Runs in your browser, nothing to install", 400, 32, MUTED, 1.05),
                  ("Genome files stay on your machine", 400, 32, MUTED, 1.2),
                  ("SVG, PNG, and PDF, plus a CLI and Python API", 400, 32, MUTED, 1.35))

    def _outro_rows(self) -> tuple[tuple[str, int, int, tuple, float], ...]:
        """The outro's rows: a story's own text in the default rows' styles, or the default rows."""

        texts = self.recording.story.outro
        if len(texts) > len(self.OUTRO_ROWS):
            raise ValueError(f"An outro has at most {len(self.OUTRO_ROWS)} rows: {texts!r}")
        return tuple((text, *row[1:]) for text, row in zip(texts, self.OUTRO_ROWS)) if texts else self.OUTRO_ROWS

    def _outro_frame(self, local: int) -> Image.Image:
        t = local / FPS
        outro_rows = self._outro_rows()
        canvas = self.background.copy().convert("RGBA")
        width, height = self.layout.size
        p = _ease(t / 0.7)
        if self.layout.portrait:
            box = self.layout.card
            card = self._outro_card((box[2], 600))
            logo = _white_logo(430)
            # Center the group on the window area, so a wide, short figure leaves no empty foot.
            rows = sum(size + (30 if weight == 700 else 16) for _, weight, size, _, _ in outro_rows)
            top = box[1] + max(0, (box[3] - (card.height + 50 + logo.height + 30 + rows)) // 2)
            canvas.alpha_composite(_fade(card, p), ((width - card.width) // 2, top + round((1 - p) * 20)))
            y = top + card.height + 50
            canvas.alpha_composite(_fade(logo, _ease((t - 0.3) / 0.7)), ((width - logo.width) // 2, y))
            y += logo.height + 30
            for text, weight, size, color, delay in outro_rows:
                self._text(canvas, self.lang(text), weight, size, color, _ease((t - delay) / 0.5), (width / 2, y),
                           center=True, rise=0)
                y += size + (30 if weight == 700 else 16)
            return canvas.convert("RGB")
        card = self._outro_card((996, 816))
        canvas.alpha_composite(_fade(card, p), (110, (height - card.height) // 2 + round((1 - p) * 20)))
        right = 110 + card.width + 90
        logo = _white_logo(560)
        canvas.alpha_composite(_fade(logo, _ease((t - 0.3) / 0.7)), (right, 250))
        y = 250 + logo.height + 50
        for text, weight, size, color, delay in outro_rows:
            self._text(canvas, self.lang(text), weight, size, color, _ease((t - delay) / 0.5), (right, y), rise=0)
            y += size + (34 if weight == 700 else 22)
        return canvas.convert("RGB")

    # timeline -------------------------------------------------------------
    def _render(self, piece: Piece, local: int, edit: Edit) -> Image.Image:
        if piece.kind == "recording":
            return self._recording_frame(piece, local, edit)
        if piece.kind == "finale":
            return self._finale_frame(piece, local, edit)
        if piece.kind == "intro":
            return self._intro_frame(local)
        if piece.kind == "montage":
            return self._montage_frame(local)
        return self._outro_frame(local)

    def frames(self, edit: Edit):
        """Yield every frame of ``edit``, dissolving across piece boundaries."""

        previous = None
        for piece in edit.pieces:
            for local in range(piece.frames):
                image = self._render(piece, local, edit)
                if previous is not None and local < piece.fade:
                    # The outgoing piece keeps playing under the dissolve.
                    under = self._render(previous, previous.frames + local, edit)
                    image = Image.blend(under, image, _ease((local + 1) / (piece.fade + 1)))
                yield image
            previous = piece

    def captions(self, edit: Edit) -> list[tuple[float, float, str]]:
        """Caption runs, in output seconds, for a subtitle file."""

        runs: list[list] = []
        position = 0
        for piece in edit.pieces:
            for local in range(0, piece.frames, 3):
                if piece.caption:
                    step, title, subtitle = piece.caption
                elif piece.kind == "recording":
                    step, title, subtitle = self.recording.caption_at(piece.moment(local), piece.start)[0]
                elif piece.kind == "finale":
                    step, title, subtitle = edit.finale_caption
                else:
                    step = title = subtitle = None
                text = None
                if title:
                    title, subtitle = self.lang(title), self.lang(subtitle)
                    head = f"{step}. {title}" if edit.numbered and step else title
                    text = head + (f"\n{subtitle}" if subtitle else "")
                second = (position + local) / FPS
                if runs and runs[-1][2] == text:
                    runs[-1][1] = second + 3 / FPS
                else:
                    if runs:
                        # A piece of a length not divisible by 3 would let the cue run past the next one's start.
                        runs[-1][1] = min(runs[-1][1], second)
                    runs.append([second, second + 3 / FPS, text])
            position += piece.frames
        return [(a, b, text) for a, b, text in runs if text]


def _srt_time(seconds: float) -> str:
    millis = round(seconds * 1000)
    return f"{millis // 3_600_000:02}:{millis // 60_000 % 60:02}:{millis // 1000 % 60:02},{millis % 1000:03}"


def _encode(frames, video: Path, review: Path, layout: Layout) -> int:
    width, height = layout.size
    command = ["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{width}x{height}", "-r", str(FPS), "-i", "-"]
    if layout.audio:
        # Some social upload paths reject or re-time video-only files.
        command += ["-f", "lavfi", "-i", "anullsrc=channel_layout=stereo:sample_rate=44100",
                    "-shortest", "-c:a", "aac", "-b:a", "128k"]
    else:
        command += ["-an"]
    command += ["-c:v", "libx264", "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p",
                "-movflags", "+faststart", str(video)]
    encoder = subprocess.Popen(command, stdin=subprocess.PIPE)
    review.mkdir(parents=True, exist_ok=True)
    count = 0
    try:
        for count, image in enumerate(frames, 1):
            encoder.stdin.write(image.tobytes())
            if (count - 1) % (FPS * 2) == 0:
                image.save(review / f"{count - 1:05d}.png")
    finally:
        encoder.stdin.close()
        encoder.wait()
    if encoder.returncode:
        raise RuntimeError(f"FFmpeg failed to encode {video.name}")
    return count


def _outputs(story: Story):
    lang = Language(story.language)
    for name, layout, stem in story.outputs:
        yield lang, name, layout, f"{stem}{lang.suffix}"


def render_walkthrough(run: Path, out: Path) -> Path:
    raw = run / "raw" / "walkthrough"
    final = out / "final"
    reports = out / "reports"
    work = out / "work"
    for path in (final, reports, work):
        path.mkdir(parents=True, exist_ok=True)
    recording = Recording(raw)
    story = recording.story
    composers: dict = {}
    report = {"recording_manifest_sha256": sha256(raw / "walkthrough.json"), "journey": recording.meta["journey"],
              "videos": {}, "covers": {}}
    for lang, name, layout, stem in _outputs(story):
        key = (lang.code, layout.name)
        if key not in composers:
            edits = {edit for edit, other, _ in story.outputs if other == layout}
            composers[key] = Composer(recording, lang.type(work), layout, lang, edits)
        composer = composers[key]
        edit = composer.edits[name]
        video = final / f"{stem}.mp4"
        frames = _encode(composer.frames(edit), video, reports / "review-frames" / stem, layout)
        if frames != sum(piece.frames for piece in edit.pieces):
            raise AssertionError(f"{stem} encoded {frames} frames")
        srt = [f"{number}\n{_srt_time(a)} --> {_srt_time(b)}\n{text}\n"
               for number, (a, b, text) in enumerate(composer.captions(edit), 1)]
        (final / f"{stem}{'.en' if lang.code == 'en' else ''}.srt").write_text("\n".join(srt), encoding="utf-8")
        timeline, position = [], 0
        for piece in edit.pieces:
            timeline.append({"kind": piece.kind, "at": round(position / FPS, 3),
                             "seconds": round(piece.frames / FPS, 3), "from": round(piece.start, 3),
                             "speed": piece.speed})
            position += piece.frames
        report["videos"][stem] = {
            "video": str(video.relative_to(out)), "video_sha256": sha256(video), "language": lang.code,
            "size": list(layout.size), "audio": layout.audio,
            "frames": frames, "seconds": round(frames / FPS, 3), "timeline": timeline,
        }
    # Every edit ends on the same outro; its last frame has faded fully in.
    last = round(OUTRO_SECONDS * FPS) - 1
    suffix = Language(story.language).suffix
    covers: list[tuple[str, Image.Image]] = []
    landscape = composers.get((story.language, "landscape"))
    if landscape is not None:
        covers.append(("poster", landscape._outro_frame(last)))
    portrait = composers.get((story.language, "portrait"))
    if portrait is not None:
        cover = portrait._outro_frame(last)
        # RedNote prefers a 3:4 cover; keep the window area and centre it.
        card = PORTRAIT.card
        top = card[1] + card[3] // 2 - 720
        covers += [("cover-vertical", cover), ("cover-3x4", cover.crop((0, top, 1080, top + 1440)))]
    for name, image in covers:
        image.save(final / f"{name}{suffix}.png")
        report["covers"][f"{name}{suffix}.png"] = list(image.size)
    report["fast_forward"] = [{"start": round(a, 2), "end": round(b, 2), "factor": round(f, 2), "kind": kind}
                              for a, b, f, kind in recording.timemap.segments]
    (reports / "walkthrough-report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return final / f"{story.outputs[0][2]}.mp4"


def check_walkthrough(out: Path) -> dict:
    """Decode the rendered videos and verify them against their report."""

    report = json.loads((out / "reports" / "walkthrough-report.json").read_text(encoding="utf-8"))
    summary = {"status": "PASS", "journey": report["journey"], "videos": {}, "fast_forward": report["fast_forward"]}
    for stem, row in report["videos"].items():
        video = out / row["video"]
        if sha256(video) != row["video_sha256"]:
            raise AssertionError(f"{video.name} does not match its render report")
        probe = json.loads(subprocess.run(
            ["ffprobe", "-v", "error", "-count_frames", "-show_entries",
             "stream=codec_type,codec_name,width,height,r_frame_rate,pix_fmt,nb_read_frames",
             "-of", "json", str(video)], check=True, capture_output=True, text=True).stdout)
        streams = {stream["codec_type"]: stream for stream in probe["streams"]}
        picture = streams.get("video", {})
        expected = {"codec_name": "h264", "width": row["size"][0], "height": row["size"][1],
                    "r_frame_rate": f"{FPS}/1", "pix_fmt": "yuv420p"}
        if any(picture.get(key) != value for key, value in expected.items()) \
                or int(picture.get("nb_read_frames", -1)) != row["frames"] \
                or set(streams) != ({"video", "audio"} if row["audio"] else {"video"}) \
                or (row["audio"] and streams["audio"].get("codec_name") != "aac"):
            raise AssertionError(f"{video.name} contract failed: {probe}")
        subprocess.run(["ffmpeg", "-v", "error", "-i", str(video), "-f", "null", "-"], check=True)
        subtitles = out / "final" / f"{stem}{'.en' if row['language'] == 'en' else ''}.srt"
        if not subtitles.is_file():
            raise AssertionError(f"Missing {subtitles.name}")
        summary["videos"][stem] = {"video": row["video"], "size": row["size"], "seconds": row["seconds"]}
    for name, size in report["covers"].items():
        path = out / "final" / name
        with Image.open(path) as image:
            if list(image.size) != size:
                raise AssertionError(f"{path.name} is {image.size}, expected {size}")
    return summary
