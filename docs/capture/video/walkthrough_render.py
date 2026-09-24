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
                if factor > 1.05:
                    segments.append((open_segment["t"], event["t"], factor))
                open_segment = None
        starts, factors, outs = [begin], [1.0], [0.0]
        for start, stop, factor in segments:
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


# -- segments ---------------------------------------------------------------

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
        self.begin = begin
        self.camera = Camera(self.events, self.timemap, self.viewport)
        self.pointer = Pointer(self.events, self.timemap)
        self.frame_times = [row[0] for row in self.frames]
        self.captions = [(self.timemap.out(e["t"]), (e["step"], e["title"], e["subtitle"]))
                         for e in self.events if e["kind"] == "caption"]
        self.toasts = [(self.timemap.out(e["t"]), (e["text"], e["icon"]))
                       for e in self.events if e["kind"] == "toast"]
        self._decoded: tuple[int, Image.Image] | None = None

    @property
    def duration(self) -> float:
        return self.timemap.duration

    def source_frame(self, output: float) -> Image.Image:
        moment = self.timemap.source(output) + self.begin
        index = max(0, bisect.bisect_right(self.frame_times, moment) - 1)
        if self._decoded is None or self._decoded[0] != index:
            with Image.open(self.raw / "frames" / self.frames[index][1]) as image:
                self._decoded = (index, image.convert("RGB"))
        return self._decoded[1]

    def caption_at(self, output: float) -> tuple[tuple, float, tuple | None]:
        """Current caption, its fade-in progress, and the caption it replaces."""

        index = bisect.bisect_right([c[0] for c in self.captions], output) - 1
        if index < 0:
            return self.captions[0][1], 0.0, None
        start, current = self.captions[index]
        previous = self.captions[index - 1][1] if index else None
        return current, (output - start) / 0.35, previous


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
        self.intro_frames = round(INTRO_SECONDS * FPS)
        self.record_frames = math.ceil(self.recording.duration * FPS)
        self.finale_frames = len(self.finale) if self.finale else 0
        self.outro_frames = round(OUTRO_SECONDS * FPS)
        self.total = self.intro_frames + self.record_frames + self.finale_frames + self.outro_frames
        self.chapters = {
            "intro": 0, "recording": self.intro_frames,
            "finale": self.intro_frames + self.record_frames,
            "outro": self.intro_frames + self.record_frames + self.finale_frames,
        }

    # frames ---------------------------------------------------------------
    def _screen(self, output: float) -> Image.Image:
        rec = self.recording
        frame = rec.source_frame(output)
        cx, cy, zoom = rec.camera.at(output)
        vw, vh = rec.viewport[0] / zoom, rec.viewport[1] / zoom
        left, top = cx - vw / 2, cy - vh / 2
        s = rec.scale
        width, height = CARD[2], CARD[3]
        screen = frame.resize((width, height), Image.Resampling.BICUBIC,
                              box=(left * s, top * s, (left + vw) * s, (top + vh) * s))
        # Pointer, drawn at a constant size after the camera transform.
        px, py = rec.pointer.at(output)
        sx, sy = (px - left) / vw * width, (py - top) / vh * height
        draw = ImageDraw.Draw(screen, "RGBA")
        for rx, ry, p in rec.pointer.ripples(output):
            ox, oy = (rx - left) / vw * width, (ry - top) / vh * height
            radius = 10 + 26 * _ease(p)
            alpha = round(150 * (1 - p))
            draw.ellipse((ox - radius, oy - radius, ox + radius, oy + radius), outline=ACCENT + (alpha,), width=4)
        cursor = self.cursor
        if rec.pointer.pressed(output):
            cursor = cursor.resize((round(cursor.width * 0.88), round(cursor.height * 0.88)), Image.Resampling.LANCZOS)
        hx = self.hotspot[0] * cursor.width / self.cursor.width
        hy = self.hotspot[1] * cursor.height / self.cursor.height
        screen.paste(cursor, (round(sx - hx), round(sy - hy)), cursor)
        return screen

    def _compose_card(self, screen: Image.Image, caption: tuple, caption_in: float,
                      previous: tuple | None, extras: list[tuple[Image.Image, tuple[int, int]]]) -> Image.Image:
        canvas = self.base.copy()
        canvas.paste(screen, CARD[:2], self.mask)
        canvas = canvas.convert("RGBA")
        for layer, position in extras:
            canvas.alpha_composite(layer, position)
        y = CARD[1] + CARD[3] + 4
        if previous is not None and caption_in < 1:
            canvas.alpha_composite(_fade(self.overlay.caption(*previous), 1 - _ease(caption_in * 2)), (CARD[0], y))
        p = _ease(caption_in)
        if p > 0:
            canvas.alpha_composite(_fade(self.overlay.caption(*caption), p), (CARD[0], y + round((1 - p) * 14)))
        brand = self.overlay.brand
        canvas.alpha_composite(brand, (CARD[0] + CARD[2] - brand.width, y + 30))
        return canvas.convert("RGB")

    def _recording_frame(self, output: float) -> Image.Image:
        rec = self.recording
        extras = []
        for start, (text, icon) in rec.toasts:
            age = output - start
            if 0 <= age < 2.8:
                alpha = min(_ease(age / 0.3), _ease((2.8 - age) / 0.4))
                layer = self.overlay.toast(text, icon)
                x = CARD[0] + (CARD[2] - layer.width) // 2
                extras.append((_fade(layer, alpha), (x, CARD[1] + CARD[3] - 140 + round((1 - alpha) * 16))))
        factor = rec.timemap.factor(output)
        if factor > 1.3:
            badge = self.overlay.badge(factor)
            extras.append((badge, (CARD[0] + CARD[2] - badge.width - 22, CARD[1] + 20)))
        caption, caption_in, previous = rec.caption_at(output)
        return self._compose_card(self._screen(output), caption, caption_in, previous, extras)

    def _finale_frame(self, index: int) -> Image.Image:
        with Image.open(self.finale[index]) as image:
            screen = image.convert("RGB").resize(CARD[2:], Image.Resampling.LANCZOS)
        progress = index / FPS
        return self._compose_card(screen, FINALE_CAPTION, progress / 0.35 + 0.0, None, [])

    def _intro_frame(self, index: int) -> Image.Image:
        t = index / FPS
        canvas = self.background.copy().convert("RGBA")
        logo = _white_logo(820)
        p = _ease(t / 0.9)
        scaled = logo.resize((round(logo.width * (0.96 + 0.04 * p)), round(logo.height * (0.96 + 0.04 * p))),
                             Image.Resampling.LANCZOS)
        canvas.alpha_composite(_fade(scaled, p), ((SIZE[0] - scaled.width) // 2, 300 - (scaled.height - logo.height) // 2))
        draw = ImageDraw.Draw(canvas)
        lines = (("From a GenBank file to a publication-ready genome map", 600, 46, INK, 0.5, 580),
                 ("A live session in the gbdraw web app. Waits are fast-forwarded and marked.", 400, 28, MUTED, 0.9, 656))
        for text, weight, size, color, delay, y in lines:
            self.type.check(text)
            font = self.type.font(weight, size)
            q = _ease((t - delay) / 0.6)
            if q <= 0:
                continue
            width = draw.textlength(text, font=font)
            layer = Image.new("RGBA", (int(width) + 4, size + 20), (0, 0, 0, 0))
            ImageDraw.Draw(layer).text((2, 0), text, font=font, fill=color + (255,))
            canvas.alpha_composite(_fade(layer, q), (round((SIZE[0] - width) / 2), y + round((1 - q) * 12)))
        return canvas.convert("RGB")

    def _outro_frame(self, index: int) -> Image.Image:
        t = index / FPS
        canvas = self.background.copy().convert("RGBA")
        figure = self.figure.copy()
        figure.thumbnail((980, 700), Image.Resampling.LANCZOS)
        card = Image.new("RGBA", (figure.width + 48, figure.height + 48), (255, 255, 255, 255))
        card.paste(figure, (24, 24))
        mask = Image.new("L", (card.width * 3, card.height * 3), 0)
        ImageDraw.Draw(mask).rounded_rectangle((0, 0, mask.width - 1, mask.height - 1), RADIUS * 3, fill=255)
        card.putalpha(mask.resize(card.size, Image.Resampling.LANCZOS))
        p = _ease(t / 0.7)
        canvas.alpha_composite(_fade(card, p), (110, (SIZE[1] - card.height) // 2 + round((1 - p) * 20)))
        right = 110 + card.width + 90
        logo = _white_logo(560)
        q = _ease((t - 0.3) / 0.7)
        canvas.alpha_composite(_fade(logo, q), (right, 250))
        draw = ImageDraw.Draw(canvas)
        rows = (("gbdraw.app", 700, 60, INK, 0.6),
                ("Free and open source (MIT)", 400, 32, MUTED, 0.9),
                ("Runs in your browser, nothing to install", 400, 32, MUTED, 1.05),
                ("Genome files stay on your machine", 400, 32, MUTED, 1.2),
                ("SVG, PNG, and PDF, plus a CLI and Python API", 400, 32, MUTED, 1.35))
        y = 250 + logo.height + 50
        for text, weight, size, color, delay in rows:
            self.type.check(text)
            font = self.type.font(weight, size)
            r = _ease((t - delay) / 0.5)
            if r > 0:
                layer = Image.new("RGBA", (900, size + 24), (0, 0, 0, 0))
                ImageDraw.Draw(layer).text((0, 0), text, font=font, fill=color + (255,))
                canvas.alpha_composite(_fade(layer, r), (right + round((1 - r) * 18), y))
            y += size + (34 if weight == 700 else 22)
        del draw
        return canvas.convert("RGB")

    def frame(self, number: int) -> Image.Image:
        chapters = self.chapters
        if number < chapters["recording"]:
            return self._intro_frame(number)
        if number < chapters["finale"]:
            return self._recording_frame((number - chapters["recording"]) / FPS)
        if number < chapters["outro"]:
            return self._finale_frame(number - chapters["finale"])
        return self._outro_frame(number - chapters["outro"])

    def frame_with_transitions(self, number: int) -> Image.Image:
        """Cross-fade across each chapter boundary."""

        image = self.frame(number)
        fade = round(FADE_SECONDS * FPS)
        for boundary in (self.chapters["recording"], self.chapters["finale"], self.chapters["outro"]):
            if boundary <= number < boundary + fade:
                previous = self.frame(boundary - 1)
                return Image.blend(previous, image, _ease((number - boundary + 1) / fade))
        return image


def _srt_time(seconds: float) -> str:
    millis = round(seconds * 1000)
    return f"{millis // 3_600_000:02}:{millis // 60_000 % 60:02}:{millis // 1000 % 60:02},{millis % 1000:03}"


def render_walkthrough(run: Path, out: Path) -> Path:
    raw = run / "raw" / "walkthrough"
    final = out / "final"
    reports = out / "reports"
    work = out / "work"
    for path in (final, reports, work):
        path.mkdir(parents=True, exist_ok=True)
    composer = Composer(raw, work)
    video = final / "gbdraw-walkthrough.mp4"
    encoder = subprocess.Popen(
        ["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
         "-s", f"{SIZE[0]}x{SIZE[1]}", "-r", str(FPS), "-i", "-", "-an", "-c:v", "libx264",
         "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(video)],
        stdin=subprocess.PIPE)
    review = reports / "review-frames"
    review.mkdir(exist_ok=True)
    picks = set(range(0, composer.total, FPS * 4))
    try:
        for number in range(composer.total):
            image = composer.frame_with_transitions(number)
            encoder.stdin.write(image.tobytes())
            if number in picks:
                image.save(review / f"{number:05d}.png")
    finally:
        encoder.stdin.close()
        encoder.wait()
    if encoder.returncode:
        raise RuntimeError("FFmpeg failed to encode the walkthrough")
    poster = composer.frame(composer.chapters["outro"] + composer.outro_frames - 1)
    poster.save(final / "poster.png")
    srt = []
    for number, (start, (step, title, subtitle)) in enumerate(composer.recording.captions, 1):
        nxt = composer.recording.captions[number][0] if number < len(composer.recording.captions) else composer.recording.duration
        a = start + composer.intro_frames / FPS
        b = nxt + composer.intro_frames / FPS
        body = f"{step}. {title}" if step else title
        srt.append(f"{number}\n{_srt_time(a)} --> {_srt_time(b)}\n{body}" + (f"\n{subtitle}" if subtitle else "") + "\n")
    (final / "gbdraw-walkthrough.en.srt").write_text("\n".join(srt), encoding="utf-8")
    report = {
        "video": str(video.relative_to(out)), "video_sha256": sha256(video),
        "frames": composer.total, "seconds": round(composer.total / FPS, 3),
        "chapters": {name: round(start / FPS, 3) for name, start in composer.chapters.items()},
        "fast_forward": [{"start": round(a, 2), "end": round(b, 2), "factor": round(f, 2)}
                         for a, b, f in composer.recording.timemap.segments],
        "recording_manifest_sha256": sha256(raw / "walkthrough.json"),
    }
    (reports / "walkthrough-report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return video


def check_walkthrough(out: Path) -> dict:
    """Decode the rendered walkthrough and verify it against its report."""

    report = json.loads((out / "reports" / "walkthrough-report.json").read_text(encoding="utf-8"))
    video = out / report["video"]
    if sha256(video) != report["video_sha256"]:
        raise AssertionError("Walkthrough MP4 does not match its render report")
    probe = json.loads(subprocess.run(
        ["ffprobe", "-v", "error", "-count_frames", "-show_entries",
         "stream=codec_type,codec_name,width,height,r_frame_rate,pix_fmt,nb_read_frames",
         "-of", "json", str(video)], check=True, capture_output=True, text=True).stdout)
    streams = probe["streams"]
    expected = {"codec_type": "video", "codec_name": "h264", "width": SIZE[0], "height": SIZE[1],
                "r_frame_rate": f"{FPS}/1", "pix_fmt": "yuv420p"}
    if len(streams) != 1 or any(streams[0].get(key) != value for key, value in expected.items()) \
            or int(streams[0]["nb_read_frames"]) != report["frames"]:
        raise AssertionError(f"Walkthrough MP4 contract failed: {probe}")
    subprocess.run(["ffmpeg", "-v", "error", "-i", str(video), "-f", "null", "-"], check=True)
    for name in ("poster.png", "gbdraw-walkthrough.en.srt"):
        if not (out / "final" / name).is_file():
            raise AssertionError(f"Missing walkthrough output: {name}")
    return {"status": "PASS", "video": report["video"], "frames": report["frames"],
            "seconds": report["seconds"], "fast_forward": report["fast_forward"]}
