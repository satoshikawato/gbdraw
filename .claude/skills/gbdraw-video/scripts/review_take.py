"""Review sheets and checks for one rendered gbdraw video.

usage: review_take.py MP4 SRT OUT_DIR [--opening SECONDS] [--min SECONDS] [--max SECONDS]
                      [--max-cue SECONDS] [--no-audio]

Writes OUT_DIR/cues.png (the start, middle and end of every caption cue, and a frame every 3 s inside a long
cue, each labeled with its cue text), OUT_DIR/opening.png with --opening (the opening every 0.25 s), and
OUT_DIR/report.json. Prints the failures and warnings, and exits 1 when a check fails.
Needs ffmpeg, ffprobe and Pillow. Cue text is drawn with NotoSansSC-Regular.otf from GBDRAW_VIDEO_FONT_DIR or
~/.cache/gbdraw-video-fonts when it is there.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

LABEL = 44  # the label band above each thumbnail
EDGE = 0.15  # how far inside a cue its first and last frames are taken


def label_font(size: int) -> ImageFont.ImageFont | ImageFont.FreeTypeFont:
    dirs = [Path(os.environ["GBDRAW_VIDEO_FONT_DIR"])] if os.environ.get("GBDRAW_VIDEO_FONT_DIR") else []
    for folder in [*dirs, Path.home() / ".cache/gbdraw-video-fonts"]:
        face = folder / "NotoSansSC-Regular.otf"
        if face.exists():
            return ImageFont.truetype(str(face), size)
    return ImageFont.load_default()


def seconds(stamp: str) -> float:
    hours, minutes, rest = stamp.replace(",", ".").split(":")
    return int(hours) * 3600 + int(minutes) * 60 + float(rest)


def read_cues(srt: Path) -> list[dict]:
    cues = []
    for block in re.split(r"\n\s*\n", srt.read_text(encoding="utf-8").strip()):
        lines = block.splitlines()
        stamp = next((line for line in lines if "-->" in line), None)
        if stamp is None:
            continue
        start, end = (seconds(part.strip()) for part in stamp.split("-->"))
        text = " / ".join(line for line in lines[lines.index(stamp) + 1:] if line.strip())
        cues.append({"index": len(cues) + 1, "start": start, "end": end, "text": text})
    return cues


def probe(mp4: Path) -> dict:
    out = subprocess.check_output(["ffprobe", "-v", "error", "-show_entries",
                                   "stream=codec_type,codec_name,width,height:format=duration",
                                   "-of", "json", str(mp4)])
    return json.loads(out)


def grab(mp4: Path, t: float, frame: Path) -> None:
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-ss", f"{t:.3f}", "-i", str(mp4), "-frames:v", "1",
                    str(frame)], check=True)


def tile(items: list[tuple[str, Path]], out: Path, size: tuple[int, int], cols: int) -> None:
    w, h = size
    rows = (len(items) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * w, rows * (h + LABEL)), "white")
    draw = ImageDraw.Draw(sheet)
    font = label_font(15)
    for i, (label, frame) in enumerate(items):
        x, y = (i % cols) * w, (i // cols) * (h + LABEL)
        sheet.paste(Image.open(frame).convert("RGB").resize((w, h)), (x, y + LABEL))
        lines, line = [], ""
        for char in label:
            if draw.textlength(line + char, font=font) > w - 8:
                lines.append(line)
                line = ""
            line += char
        lines.append(line)
        for n, text in enumerate(lines[:2]):
            draw.text((x + 4, y + 3 + n * 19), text, fill="black", font=font)
    sheet.save(out)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mp4", type=Path)
    parser.add_argument("srt", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--opening", type=float, default=0.0, help="seconds of opening to sample every 0.25 s")
    parser.add_argument("--min", type=float, help="shortest allowed length in seconds")
    parser.add_argument("--max", type=float, help="longest allowed length in seconds")
    parser.add_argument("--max-cue", type=float, default=12.0, help="warn about cues longer than this")
    parser.add_argument("--no-audio", action="store_true", help="fail when the video has an audio stream")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    frames = args.out / "frames"
    frames.mkdir(exist_ok=True)
    failures: list[str] = []
    warnings: list[str] = []

    info = probe(args.mp4)
    duration = float(info["format"]["duration"])
    video = [s for s in info["streams"] if s["codec_type"] == "video"]
    if len(video) != 1:
        failures.append(f"{len(video)} video streams")
    size = (video[0]["width"], video[0]["height"]) if video else (1080, 1920)
    if size not in ((1080, 1920), (1920, 1080)):
        failures.append(f"size {size[0]}x{size[1]}, not 1080x1920 or 1920x1080")
    if video and video[0]["codec_name"] != "h264":
        failures.append(f"codec {video[0]['codec_name']}, not h264")
    if args.no_audio and any(s["codec_type"] == "audio" for s in info["streams"]):
        failures.append("an audio stream")
    if args.min is not None and duration < args.min:
        failures.append(f"length {duration:.1f} s < {args.min} s")
    if args.max is not None and duration > args.max:
        failures.append(f"length {duration:.1f} s > {args.max} s")

    cues = read_cues(args.srt)
    for before, cue in zip(cues, cues[1:]):
        if cue["start"] < before["end"] - 0.001:
            failures.append(f"cue {cue['index']} starts before cue {before['index']} ends")
    for cue in cues:
        if cue["end"] > duration + 0.05:
            failures.append(f"cue {cue['index']} ends after the video")
        if cue["end"] - cue["start"] > args.max_cue:
            warnings.append(f"cue {cue['index']} lasts {cue['end'] - cue['start']:.1f} s: check that its text is "
                            f"true for every frame ({cue['text']})")

    portrait = size[0] < size[1]
    items = []
    for cue in cues:
        start, end = cue["start"] + EDGE, max(cue["start"] + EDGE, cue["end"] - EDGE)
        times = sorted({start, (cue["start"] + cue["end"]) / 2, end,
                        *(t for t in (cue["start"] + 3.0 * k for k in range(1, 100)) if t < end - 1.0)})
        for t in times:
            frame = frames / f"cue{cue['index']:02d}-{t:07.2f}.png"
            grab(args.mp4, t, frame)
            items.append((f"{t:.1f}s #{cue['index']} {cue['text']}", frame))
    if items:
        tile(items, args.out / "cues.png", (270, 480) if portrait else (480, 270), 8 if portrait else 4)
    if args.opening > 0:
        shots = []
        t = 0.125
        while t < min(args.opening, duration):
            frame = frames / f"opening-{t:06.3f}.png"
            grab(args.mp4, t, frame)
            shots.append((f"{t:.2f}s", frame))
            t += 0.25
        tile(shots, args.out / "opening.png", (360, 640) if portrait else (640, 360), 8 if portrait else 4)

    report = {"video": str(args.mp4), "duration": round(duration, 3), "size": list(size), "cues": len(cues),
              "failures": failures, "warnings": warnings}
    (args.out / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{args.mp4.name}: {duration:.1f} s, {size[0]}x{size[1]}, {len(cues)} cues; sheets in {args.out}")
    for line in failures:
        print("FAIL", line)
    for line in warnings:
        print("WARN", line)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
