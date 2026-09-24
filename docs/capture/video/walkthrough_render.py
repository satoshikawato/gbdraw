"""Compose a recorded walkthrough into a captioned, camera-guided MP4.

Inputs are the run-local files written by ``walkthrough.record_walkthrough``.
This module never starts a browser, so captions, pacing, and camera cues can be
revised from one recording.
"""

from __future__ import annotations

import bisect
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

from fontTools.ttLib import TTFont
from PIL import Image, ImageDraw, ImageFilter, ImageFont

from video.model import sha256

REPO_ROOT = Path(__file__).resolve().parents[3]
INTER = REPO_ROOT / "gbdraw/web/vendor/fonts/inter"
LOGO = REPO_ROOT / "gbdraw/web/assets/gbdraw-logo-title.png"

FPS = 30
SIZE = (1920, 1080)
CARD = (160, 34, 1600, 900)  # x, y, width, height of the app window
RADIUS = 18
INTRO_SECONDS = 3.4
OUTRO_SECONDS = 5.0
FADE_SECONDS = 0.5
FINALE_CAPTION = (7, "Vector output stays sharp at any zoom", "The downloaded SVG, magnified 7×")
INK = (241, 245, 249)
MUTED = (148, 170, 196)
ACCENT = (59, 130, 246)


def _ease(value: float) -> float:
    value = min(1.0, max(0.0, value))
    return value * value * (3 - 2 * value)


def _fonts(work: Path) -> dict[str, Path]:
    """Convert the app's vendored Inter WOFF2 files once per render."""

    fonts = {}
    for weight in (400, 500, 600, 700):
        target = work / f"inter-{weight}.ttf"
        if not target.exists():
            font = TTFont(INTER / f"inter-latin-{weight}-normal.woff2")
            font.flavor = None
            font.save(target)
        fonts[weight] = target
    return fonts


class Type:
    def __init__(self, fonts: dict[str, Path]) -> None:
        self._fonts = fonts
        self._cache: dict = {}
        self._cmap = set(TTFont(fonts[400]).getBestCmap())

    def font(self, weight: int, size: int) -> ImageFont.FreeTypeFont:
        key = (weight, size)
        if key not in self._cache:
            self._cache[key] = ImageFont.truetype(str(self._fonts[weight]), size)
        return self._cache[key]

    def check(self, text: str) -> None:
        missing = {char for char in text if ord(char) not in self._cmap and char not in "\n"}
        if missing:
            raise ValueError(f"Inter latin subset lacks glyphs {sorted(missing)!r} in {text!r}")


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
    """Ease between logged framing targets, clamped to the page."""

    def __init__(self, events: list[dict], timemap: TimeMap, viewport: tuple[int, int]) -> None:
        self.viewport = viewport
        self.cues = [(timemap.out(e["t"]), self._clamp(e["cx"], e["cy"], e["zoom"]), e["duration"])
                     for e in events if e["kind"] == "camera"]
        self._states = []
        state = self._clamp(viewport[0] / 2, viewport[1] / 2, 1.0)
        for index, (start, target, duration) in enumerate(self.cues):
            if index:
                state = self._at(index - 1, start)
            self._states.append(state)

    def _clamp(self, cx: float, cy: float, zoom: float) -> tuple[float, float, float]:
        zoom = max(1.0, zoom)
        half_w, half_h = self.viewport[0] / zoom / 2, self.viewport[1] / zoom / 2
        return (min(max(cx, half_w), self.viewport[0] - half_w),
                min(max(cy, half_h), self.viewport[1] - half_h), zoom)

    def _at(self, index: int, output: float) -> tuple[float, float, float]:
        start, target, duration = self.cues[index]
        origin = self._states[index]
        p = 1.0 if duration <= 0 else _ease((output - start) / duration)
        # Interpolate zoom geometrically so zooming feels even.
        zoom = math.exp(math.log(origin[2]) + (math.log(target[2]) - math.log(origin[2])) * p)
        return (origin[0] + (target[0] - origin[0]) * p, origin[1] + (target[1] - origin[1]) * p, zoom)

    def at(self, output: float) -> tuple[float, float, float]:
        starts = [cue[0] for cue in self.cues]
        index = bisect.bisect_right(starts, output) - 1
        if index < 0:
            return self._clamp(self.viewport[0] / 2, self.viewport[1] / 2, 1.0)
        return self._clamp(*self._at(index, output))


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

def _background() -> Image.Image:
    width, height = SIZE
    top, bottom = (15, 23, 42), (23, 43, 77)
    column = Image.new("RGB", (1, height))
    for y in range(height):
        p = y / (height - 1)
        column.putpixel((0, y), tuple(round(a + (b - a) * p) for a, b in zip(top, bottom)))
    image = column.resize(SIZE)
    glow = Image.new("L", SIZE, 0)
    ImageDraw.Draw(glow).ellipse((-400, -700, 1400, 700), fill=70)
    glow = glow.filter(ImageFilter.GaussianBlur(220))
    return Image.composite(Image.new("RGB", SIZE, (37, 99, 235)), image, glow.point(lambda v: v // 3))


def _card_layers(background: Image.Image) -> tuple[Image.Image, Image.Image]:
    x, y, width, height = CARD
    shadow = Image.new("L", SIZE, 0)
    ImageDraw.Draw(shadow).rounded_rectangle((x, y + 14, x + width, y + height + 14), RADIUS, fill=150)
    shadow = shadow.filter(ImageFilter.GaussianBlur(26))
    base = Image.composite(Image.new("RGB", SIZE, (2, 6, 23)), background, shadow)
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


def _white_logo(width: int) -> Image.Image:
    logo = Image.open(LOGO).convert("RGBA")
    alpha = logo.getchannel("A")
    white = Image.new("RGBA", logo.size, (255, 255, 255, 0))
    white.putalpha(alpha)
    height = round(logo.height * width / logo.width)
    return white.resize((width, height), Image.Resampling.LANCZOS)


class Overlay:
    """Pre-rendered captions, toasts, and badges, blended per frame."""

    def __init__(self, type_: Type) -> None:
        self.type = type_
        self._captions: dict = {}
        self._toasts: dict = {}
        self._badges: dict = {}
        self.brand = self._brand()

    def _brand(self) -> Image.Image:
        logo = _white_logo(150)
        font = self.type.font(500, 22)
        text = "gbdraw.app"
        width = logo.width
        layer = Image.new("RGBA", (width + 4, logo.height + 34), (0, 0, 0, 0))
        layer.alpha_composite(logo, (0, 0))
        draw = ImageDraw.Draw(layer)
        text_width = draw.textlength(text, font=font)
        draw.text(((width - text_width) / 2 + 2, logo.height + 6), text, font=font, fill=MUTED + (255,))
        return layer

    def caption(self, step: int | None, title: str, subtitle: str) -> Image.Image:
        key = (step, title, subtitle)
        if key in self._captions:
            return self._captions[key]
        self.type.check(title + subtitle)
        layer = Image.new("RGBA", (1440, 144), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        x = 0
        if step is not None:
            draw.ellipse((0, 24, 56, 80), fill=ACCENT + (255,))
            number = self.type.font(700, 30)
            label = str(step)
            width = draw.textlength(label, font=number)
            draw.text((28 - width / 2, 32), label, font=number, fill=(255, 255, 255, 255))
            x = 76
        draw.text((x, 18 if subtitle else 30), title, font=self.type.font(600, 40), fill=INK + (255,))
        if subtitle:
            draw.text((x, 74), subtitle, font=self.type.font(400, 26), fill=MUTED + (255,))
        self._captions[key] = layer
        return layer

    def toast(self, text: str, icon: str) -> Image.Image:
        key = (text, icon)
        if key in self._toasts:
            return self._toasts[key]
        self.type.check(text)
        font = self.type.font(500, 26)
        probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
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
        ImageDraw.Draw(layer).text((90, 33), text, font=font, fill=(15, 23, 42, 255))
        self._toasts[key] = layer
        return layer

    def badge(self, factor: float) -> Image.Image:
        label = f"▶▶  {factor:.0f}×" if factor >= 2 else "▶▶  fast"
        if label in self._badges:
            return self._badges[label]
        font = self.type.font(600, 24)
        text = label.replace("▶▶  ", "")
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
        self._badges[label] = layer
        return layer


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
        draw.text((left + 12 * scale, y + row / 2 - 1 * scale), label, font=font, anchor="lm",
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
        self.frames = json.loads((raw / "frames.json").read_text(encoding="utf-8"))
        self.events = json.loads((raw / "events.json").read_text(encoding="utf-8"))
        self.viewport = tuple(self.meta["viewport"])
        self.scale = self.meta["scale"]
        begin = self.frames[0][0]
        end = next(e["t"] for e in self.events if e["kind"] == "end")
        self.timemap = TimeMap.build(self.events, begin, end)
        self.camera = Camera(self.events, self.timemap, self.viewport)
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

    def source_frame(self, moment: float) -> Image.Image:
        index = max(0, bisect.bisect_right(self.frame_times, self.timemap.source(moment)) - 1)
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
    ("upload", -0.3, "loaded", 1.2),
    ("generate", 0.0, "generated", 2.0),
    ("labels", -0.2, "crowded", 2.2),
    ("priority-type", -0.2, "genes", 2.6),
    ("rule-4", -0.1, "legend", 2.4),
    ("bracket", -0.2, "dloop", 2.6),
    ("export", -0.1, "downloaded", 1.4),
)
HIGHLIGHT_SPEED = 1.4
HIGHLIGHT_WAIT_SPEED = 2.5  # generation waits, on top of the walkthrough's own compression
MONTAGE_SECONDS = 3.6
HIGHLIGHT_OUTRO_SECONDS = 4.2
CUT_FADE = 7
CHAPTER_FADE = round(FADE_SECONDS * FPS)


def full_edit(recording: Recording, finale_frames: int) -> Edit:
    return Edit("walkthrough", [
        Piece("intro", round(INTRO_SECONDS * FPS)),
        Piece("recording", math.ceil(recording.duration * FPS),
              moments=tuple(n / FPS for n in range(math.ceil(recording.duration * FPS))), fade=CHAPTER_FADE),
        Piece("finale", finale_frames, fade=CHAPTER_FADE),
        Piece("outro", round(OUTRO_SECONDS * FPS), fade=CHAPTER_FADE),
    ], numbered=True, badge_from=1.3, finale_caption=FINALE_CAPTION)


def highlights_edit(recording: Recording, finale_frames: int) -> Edit:
    pieces = [Piece("montage", round(MONTAGE_SECONDS * FPS))]
    for index, (first, before, last, after) in enumerate(HIGHLIGHTS):
        missing = {first, last} - set(recording.marks)
        if missing:
            raise ValueError(f"Recording lacks marks {sorted(missing)}")
        start = recording.marks[first] + before
        stop = recording.marks[last] + after
        if stop <= start:
            raise ValueError(f"Highlight {first}->{last} is empty")
        moments = _moments(recording.timemap, start, stop, HIGHLIGHT_SPEED, HIGHLIGHT_WAIT_SPEED)
        pieces.append(Piece("recording", len(moments), moments=moments, speed=HIGHLIGHT_SPEED,
                            fade=CHAPTER_FADE if index == 0 else CUT_FADE))
    pieces.append(Piece("finale", round(finale_frames / 1.2), speed=1.2, fade=CHAPTER_FADE))
    pieces.append(Piece("outro", round(HIGHLIGHT_OUTRO_SECONDS * FPS), fade=CHAPTER_FADE))
    return Edit("highlights", pieces, numbered=False, badge_from=2.0,
                finale_caption=(None, *FINALE_CAPTION[1:]))


# -- composition ------------------------------------------------------------

class Composer:
    def __init__(self, raw: Path, work: Path) -> None:
        self.recording = Recording(raw)
        self.type = Type(_fonts(work))
        self.overlay = Overlay(self.type)
        self.background = _background()
        self.base, self.mask = _card_layers(self.background)
        self.cursor, self.hotspot = _cursor()
        self.finale = sorted((raw / "finale").glob("*.png"))
        self.figure = Image.open(raw / "final-figure.png").convert("RGB")
        self.gallery = [(Image.open(raw / row["path"]).convert("RGB"), row["label"])
                        for row in self.recording.meta.get("gallery", [])]
        self.edits = {
            "walkthrough": full_edit(self.recording, len(self.finale)),
            "highlights": highlights_edit(self.recording, len(self.finale)),
        }
        self._static: dict = {}

    # recording ------------------------------------------------------------
    def _screen(self, moment: float) -> Image.Image:
        rec = self.recording
        frame = rec.source_frame(moment)
        cx, cy, zoom = rec.camera.at(moment)
        vw, vh = rec.viewport[0] / zoom, rec.viewport[1] / zoom
        left, top = cx - vw / 2, cy - vh / 2
        s = rec.scale
        width, height = CARD[2], CARD[3]
        screen = frame.resize((width, height), Image.Resampling.BICUBIC,
                              box=(left * s, top * s, (left + vw) * s, (top + vh) * s))
        scale = width / vw

        def to_screen(x: float, y: float) -> tuple[float, float]:
            return (x - left) * scale, (y - top) * scale

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

    def _progress(self, canvas: Image.Image, step: int | None) -> None:
        steps = self.recording.steps
        if step is None or not steps:
            return
        draw = ImageDraw.Draw(canvas)
        gap, y = 8, SIZE[1] - 14
        width = (CARD[2] - gap * (len(steps) - 1)) / len(steps)
        for index, value in enumerate(steps):
            x = CARD[0] + index * (width + gap)
            color = ACCENT if value == step else (INK if value < step else (51, 65, 85))
            draw.rounded_rectangle((x, y, x + width, y + 4), 2, fill=color)

    def _compose_card(self, screen: Image.Image, caption: tuple, caption_in: float, previous: tuple | None,
                      extras: list, edit: Edit) -> Image.Image:
        canvas = self.base.copy()
        canvas.paste(screen, CARD[:2], self.mask)
        canvas = canvas.convert("RGBA")
        for layer, position in extras:
            canvas.alpha_composite(layer, position)
        if not edit.numbered:
            caption = (None, *caption[1:])
            previous = previous and (None, *previous[1:])
        y = CARD[1] + CARD[3] + 4
        if previous is not None and caption_in < 1:
            canvas.alpha_composite(_fade(self.overlay.caption(*previous), 1 - _ease(caption_in * 2)), (CARD[0], y))
        p = _ease(caption_in)
        if p > 0:
            canvas.alpha_composite(_fade(self.overlay.caption(*caption), p), (CARD[0], y + round((1 - p) * 14)))
        brand = self.overlay.brand
        canvas.alpha_composite(brand, (CARD[0] + CARD[2] - brand.width, y + 30))
        if edit.numbered:
            self._progress(canvas, caption[0])
        return canvas.convert("RGB")

    def _recording_frame(self, piece: Piece, local: int, edit: Edit) -> Image.Image:
        rec = self.recording
        moment = piece.moment(local)
        extras = []
        for start, (text, icon) in rec.toasts:
            age = moment - start
            if 0 <= age < 2.8:
                alpha = min(_ease(age / 0.3), _ease((2.8 - age) / 0.4))
                layer = self.overlay.toast(text, icon)
                x = CARD[0] + (CARD[2] - layer.width) // 2
                extras.append((_fade(layer, alpha), (x, CARD[1] + CARD[3] - 140 + round((1 - alpha) * 16))))
        factor = rec.timemap.factor(moment) * piece.rate(local)
        if factor >= edit.badge_from:
            badge = self.overlay.badge(factor)
            extras.append((badge, (CARD[0] + CARD[2] - badge.width - 22, CARD[1] + 20)))
        caption, caption_in, previous = rec.caption_at(moment, piece.start)
        return self._compose_card(self._screen(moment), caption, caption_in, previous, extras, edit)

    def _finale_frame(self, piece: Piece, local: int, edit: Edit) -> Image.Image:
        index = min(len(self.finale) - 1, round(local * piece.speed))
        with Image.open(self.finale[index]) as image:
            screen = image.convert("RGB").resize(CARD[2:], Image.Resampling.LANCZOS)
        return self._compose_card(screen, edit.finale_caption, local / FPS / 0.35, None, [], edit)

    # cards ----------------------------------------------------------------
    def _text(self, canvas: Image.Image, text: str, weight: int, size: int, color: tuple, alpha: float,
              position: tuple[float, float], *, center: bool = False, rise: int = 12) -> None:
        if alpha <= 0:
            return
        self.type.check(text)
        font = self.type.font(weight, size)
        key = ("text", text, weight, size, color)
        if key not in self._static:
            width = ImageDraw.Draw(Image.new("RGBA", (1, 1))).textlength(text, font=font)
            layer = Image.new("RGBA", (int(width) + 4, size + 24), (0, 0, 0, 0))
            ImageDraw.Draw(layer).text((2, 0), text, font=font, fill=color + (255,))
            self._static[key] = layer
        layer = self._static[key]
        x, y = position
        if center:
            x -= layer.width / 2
        canvas.alpha_composite(_fade(layer, alpha), (round(x), round(y + (1 - alpha) * rise)))

    def _intro_frame(self, local: int) -> Image.Image:
        t = local / FPS
        canvas = self.background.copy().convert("RGBA")
        logo = _white_logo(820)
        p = _ease(t / 0.9)
        scaled = logo.resize((round(logo.width * (0.96 + 0.04 * p)), round(logo.height * (0.96 + 0.04 * p))),
                             Image.Resampling.LANCZOS)
        canvas.alpha_composite(_fade(scaled, p), ((SIZE[0] - scaled.width) // 2, 300 - (scaled.height - logo.height) // 2))
        self._text(canvas, "From a GenBank file to a publication-ready genome map", 600, 46, INK,
                   _ease((t - 0.5) / 0.6), (SIZE[0] / 2, 580), center=True)
        self._text(canvas, "A live session in the gbdraw web app. Waits are fast-forwarded and marked.", 400, 28,
                   MUTED, _ease((t - 0.9) / 0.6), (SIZE[0] / 2, 656), center=True)
        return canvas.convert("RGB")

    def _gallery_card(self, index: int, size: tuple[int, int]) -> Image.Image:
        key = ("gallery", index, size)
        if key not in self._static:
            figure, label = self.gallery[index]
            width, height = size
            card = Image.new("RGBA", size, (255, 255, 255, 255))
            fitted = figure.copy()
            fitted.thumbnail((width - 36, height - 76), Image.Resampling.LANCZOS)
            card.paste(fitted, ((width - fitted.width) // 2, 18 + (height - 76 - fitted.height) // 2))
            draw = ImageDraw.Draw(card)
            self.type.check(label)
            draw.text((22, height - 46), label, font=self.type.font(600, 24), fill=(15, 23, 42, 255))
            mask = Image.new("L", (width * 3, height * 3), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, width * 3 - 1, height * 3 - 1), 14 * 3, fill=255)
            card.putalpha(mask.resize(size, Image.Resampling.LANCZOS))
            self._static[key] = card
        return self._static[key]

    def _montage_frame(self, local: int) -> Image.Image:
        t = local / FPS
        canvas = self.background.copy().convert("RGBA")
        columns, rows = 3, 2
        gap = 28
        width = (CARD[2] - gap * (columns - 1)) // columns
        height = (CARD[3] - gap * (rows - 1)) // rows
        for index in range(min(len(self.gallery), columns * rows)):
            p = _ease((t + 0.25 - index * 0.14) / 0.45)
            if p <= 0:
                continue
            card = self._gallery_card(index, (width, height))
            grow = 0.94 + 0.06 * p
            scaled = card.resize((round(width * grow), round(height * grow)), Image.Resampling.BILINEAR)
            x = CARD[0] + (index % columns) * (width + gap) + (width - scaled.width) // 2
            y = CARD[1] + (index // columns) * (height + gap) + (height - scaled.height) // 2
            canvas.alpha_composite(_fade(scaled, p), (x, y))
        y = CARD[1] + CARD[3] + 4
        caption = self.overlay.caption(None, "Genome diagrams for microbes and organelles",
                                       "Highlights from a live session in the gbdraw web app, sped up")
        canvas.alpha_composite(_fade(caption, _ease((t - 0.8) / 0.5)), (CARD[0], y))
        canvas.alpha_composite(self.overlay.brand, (CARD[0] + CARD[2] - self.overlay.brand.width, y + 30))
        return canvas.convert("RGB")

    def _outro_frame(self, local: int) -> Image.Image:
        t = local / FPS
        canvas = self.background.copy().convert("RGBA")
        key = ("outro-card",)
        if key not in self._static:
            figure = self.figure.copy()
            figure.thumbnail((940, 760), Image.Resampling.LANCZOS)
            card = Image.new("RGBA", (figure.width + 56, figure.height + 56), (255, 255, 255, 255))
            card.paste(figure, (28, 28))
            mask = Image.new("L", (card.width * 3, card.height * 3), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, mask.width - 1, mask.height - 1), RADIUS * 3, fill=255)
            card.putalpha(mask.resize(card.size, Image.Resampling.LANCZOS))
            self._static[key] = card
        card = self._static[key]
        p = _ease(t / 0.7)
        canvas.alpha_composite(_fade(card, p), (110, (SIZE[1] - card.height) // 2 + round((1 - p) * 20)))
        right = 110 + card.width + 90
        logo = _white_logo(560)
        canvas.alpha_composite(_fade(logo, _ease((t - 0.3) / 0.7)), (right, 250))
        rows = (("gbdraw.app", 700, 60, INK, 0.6),
                ("Free and open source (MIT)", 400, 32, MUTED, 0.9),
                ("Runs in your browser, nothing to install", 400, 32, MUTED, 1.05),
                ("Genome files stay on your machine", 400, 32, MUTED, 1.2),
                ("SVG, PNG, and PDF, plus a CLI and Python API", 400, 32, MUTED, 1.35))
        y = 250 + logo.height + 50
        for text, weight, size, color, delay in rows:
            self._text(canvas, text, weight, size, color, _ease((t - delay) / 0.5), (right, y), rise=0)
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
                if piece.kind == "recording":
                    step, title, subtitle = self.recording.caption_at(piece.moment(local), piece.start)[0]
                elif piece.kind == "finale":
                    step, title, subtitle = edit.finale_caption
                else:
                    step = title = subtitle = None
                text = None
                if title:
                    head = f"{step}. {title}" if edit.numbered and step else title
                    text = head + (f"\n{subtitle}" if subtitle else "")
                second = (position + local) / FPS
                if runs and runs[-1][2] == text:
                    runs[-1][1] = second + 3 / FPS
                else:
                    runs.append([second, second + 3 / FPS, text])
            position += piece.frames
        return [(a, b, text) for a, b, text in runs if text]


def _srt_time(seconds: float) -> str:
    millis = round(seconds * 1000)
    return f"{millis // 3_600_000:02}:{millis // 60_000 % 60:02}:{millis // 1000 % 60:02},{millis % 1000:03}"


def _encode(frames, video: Path, review: Path) -> int:
    encoder = subprocess.Popen(
        ["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
         "-s", f"{SIZE[0]}x{SIZE[1]}", "-r", str(FPS), "-i", "-", "-an", "-c:v", "libx264",
         "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(video)],
        stdin=subprocess.PIPE)
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


OUTPUTS = {"walkthrough": "gbdraw-walkthrough", "highlights": "gbdraw-highlights"}


def render_walkthrough(run: Path, out: Path) -> Path:
    raw = run / "raw" / "walkthrough"
    final = out / "final"
    reports = out / "reports"
    work = out / "work"
    for path in (final, reports, work):
        path.mkdir(parents=True, exist_ok=True)
    composer = Composer(raw, work)
    report = {"recording_manifest_sha256": sha256(raw / "walkthrough.json"), "videos": {}}
    for name, stem in OUTPUTS.items():
        edit = composer.edits[name]
        video = final / f"{stem}.mp4"
        frames = _encode(composer.frames(edit), video, reports / "review-frames" / name)
        if frames != sum(piece.frames for piece in edit.pieces):
            raise AssertionError(f"{name} encoded {frames} frames")
        srt = [f"{number}\n{_srt_time(a)} --> {_srt_time(b)}\n{text}\n"
               for number, (a, b, text) in enumerate(composer.captions(edit), 1)]
        (final / f"{stem}.en.srt").write_text("\n".join(srt), encoding="utf-8")
        timeline, position = [], 0
        for piece in edit.pieces:
            timeline.append({"kind": piece.kind, "at": round(position / FPS, 3),
                             "seconds": round(piece.frames / FPS, 3), "from": round(piece.start, 3),
                             "speed": piece.speed})
            position += piece.frames
        report["videos"][name] = {
            "video": str(video.relative_to(out)), "video_sha256": sha256(video),
            "frames": frames, "seconds": round(frames / FPS, 3), "timeline": timeline,
        }
    walkthrough = composer.edits["walkthrough"]
    last = walkthrough.pieces[-1]
    composer._outro_frame(last.frames - 1).save(final / "poster.png")
    report["fast_forward"] = [{"start": round(a, 2), "end": round(b, 2), "factor": round(f, 2), "kind": kind}
                              for a, b, f, kind in composer.recording.timemap.segments]
    (reports / "walkthrough-report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return final / f"{OUTPUTS['walkthrough']}.mp4"


def check_walkthrough(out: Path) -> dict:
    """Decode the rendered videos and verify them against their report."""

    report = json.loads((out / "reports" / "walkthrough-report.json").read_text(encoding="utf-8"))
    summary = {"status": "PASS", "videos": {}, "fast_forward": report["fast_forward"]}
    for name, stem in OUTPUTS.items():
        row = report["videos"][name]
        video = out / row["video"]
        if sha256(video) != row["video_sha256"]:
            raise AssertionError(f"{video.name} does not match its render report")
        probe = json.loads(subprocess.run(
            ["ffprobe", "-v", "error", "-count_frames", "-show_entries",
             "stream=codec_type,codec_name,width,height,r_frame_rate,pix_fmt,nb_read_frames",
             "-of", "json", str(video)], check=True, capture_output=True, text=True).stdout)
        streams = probe["streams"]
        expected = {"codec_type": "video", "codec_name": "h264", "width": SIZE[0], "height": SIZE[1],
                    "r_frame_rate": f"{FPS}/1", "pix_fmt": "yuv420p"}
        if len(streams) != 1 or any(streams[0].get(key) != value for key, value in expected.items()) \
                or int(streams[0]["nb_read_frames"]) != row["frames"]:
            raise AssertionError(f"{video.name} contract failed: {probe}")
        subprocess.run(["ffmpeg", "-v", "error", "-i", str(video), "-f", "null", "-"], check=True)
        if not (out / "final" / f"{stem}.en.srt").is_file():
            raise AssertionError(f"Missing {stem}.en.srt")
        summary["videos"][name] = {"video": row["video"], "seconds": row["seconds"]}
    if not (out / "final" / "poster.png").is_file():
        raise AssertionError("Missing poster.png")
    return summary
