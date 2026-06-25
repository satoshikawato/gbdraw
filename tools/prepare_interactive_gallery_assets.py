#!/usr/bin/env python3
from __future__ import annotations

import io
import json
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cairosvg
from PIL import Image, ImageDraw, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
EXAMPLE_ROOT = GALLERY_ROOT / "examples"
SESSION_ROOT = GALLERY_ROOT / "sessions"
SOURCE_ROOT = GALLERY_ROOT / "sources"
THUMBNAIL_ROOT = GALLERY_ROOT / "thumbnails"

GENOME_SUFFIXES = (".gb", ".gbk", ".gbff")


@dataclass(frozen=True)
class GallerySessionExample:
    id: str
    title: str
    tags: tuple[str, ...]
    description: str
    interactive_step: str = "Loaded from the saved interactive SVG output paired with this session JSON."
    source_note: str = "Session JSON and generated SVG output are stored with the gallery assets."

    @property
    def session_path(self) -> Path:
        return SESSION_ROOT / f"{self.id}.gbdraw-session.json"

    @property
    def session_ref(self) -> str:
        return f"./sessions/{self.id}.gbdraw-session.json"

    @property
    def source_svg_path(self) -> Path:
        return SOURCE_ROOT / f"{self.id}.svg"

    @property
    def output_svg_path(self) -> Path:
        return self.gallery_svg_path

    @property
    def gallery_svg_path(self) -> Path:
        return EXAMPLE_ROOT / f"{self.id}.svg"

    @property
    def gallery_svg_ref(self) -> str:
        return f"./examples/{self.id}.svg"

    @property
    def thumbnail_path(self) -> Path:
        return THUMBNAIL_ROOT / f"{self.id}.webp"

    @property
    def thumbnail_ref(self) -> str:
        return f"./thumbnails/{self.id}.webp"


EXAMPLES: tuple[GallerySessionExample, ...] = (
    GallerySessionExample(
        id="Vnig_TUMSAT-TG-2018",
        title="Vibrio nigripulchritudo TUMSAT-TG-2018 complete genome",
        tags=("Circular", "Multi-record", "Session"),
        description="A multi-replicon Vibrio nigripulchritudo TUMSAT-TG-2018 circular genome diagram generated from a saved session.",
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_collinear",
        title="Hepatoplasmataceae collinear session output",
        tags=("Linear", "LOSAT", "Collinear"),
        description="A five-genome Hepatoplasmataceae comparison generated from a saved collinear BLASTP session.",
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_orthogroup",
        title="Hepatoplasmataceae orthogroup session output",
        tags=("Linear", "LOSAT", "Orthogroup"),
        description="A five-genome Hepatoplasmataceae comparison generated from a saved orthogroup BLASTP session.",
    ),
    GallerySessionExample(
        id="majanivirus_orthogroup",
        title="Majanivirus orthogroup session output",
        tags=("Linear", "LOSAT", "Orthogroup", "Large"),
        description="A nine-genome majanivirus comparison generated from a saved orthogroup BLASTP session.",
    ),
)


def _format_size(num_bytes: int) -> str:
    if num_bytes >= 1024 * 1024:
        return f"{num_bytes / (1024 * 1024):.1f} MB"
    return f"{num_bytes / 1024:.0f} KB"


def _load_session(example: GallerySessionExample) -> dict[str, Any]:
    with example.session_path.open(encoding="utf-8") as handle:
        session = json.load(handle)
    if not isinstance(session, dict):
        raise ValueError(f"{example.session_path} did not contain a JSON object.")
    return session


def _session_cli_invocation(session: dict[str, Any]) -> dict[str, Any]:
    cli = session.get("cliInvocation")
    if isinstance(cli, dict):
        return cli
    config = session.get("config")
    if isinstance(config, dict):
        cli = config.get("cliInvocation")
        if isinstance(cli, dict):
            return cli
    return {}


def _session_raw_args(session: dict[str, Any]) -> list[str]:
    cli = _session_cli_invocation(session)
    args = cli.get("args")
    if isinstance(args, list):
        return [str(arg) for arg in args]

    config = session.get("config")
    if isinstance(config, dict):
        cli_options = config.get("cliOptions")
        if isinstance(cli_options, dict):
            raw_args = cli_options.get("rawArgs")
            if isinstance(raw_args, list):
                return [str(arg) for arg in raw_args]
    return []


def _session_command(session: dict[str, Any]) -> str:
    cli = _session_cli_invocation(session)
    mode = str(cli.get("mode") or "linear")
    return shlex.join(["gbdraw", mode, *_session_raw_args(session)])


def _session_feature_sources(session: dict[str, Any]) -> list[str]:
    cli = _session_cli_invocation(session)
    seen: set[str] = set()
    sources: list[str] = []

    bindings = cli.get("fileBindings")
    if isinstance(bindings, list):
        for binding in bindings:
            if not isinstance(binding, dict):
                continue
            name = binding.get("name")
            if not isinstance(name, str) or not name.lower().endswith(GENOME_SUFFIXES):
                continue
            if name not in seen:
                seen.add(name)
                sources.append(name)
    return sources


def _remove_stale_assets() -> None:
    expected_svgs = {f"{example.id}.svg" for example in EXAMPLES}
    expected_thumbnails = {f"{example.id}.webp" for example in EXAMPLES}

    for path in EXAMPLE_ROOT.glob("*.svg"):
        if path.name not in expected_svgs:
            path.unlink()
    for path in THUMBNAIL_ROOT.glob("*.webp"):
        if path.name not in expected_thumbnails:
            path.unlink()


def _render_thumbnail(example: GallerySessionExample) -> None:
    source_path = example.source_svg_path if example.source_svg_path.exists() else example.output_svg_path
    try:
        png_bytes = cairosvg.svg2png(
            url=str(source_path),
            output_width=720,
            background_color="white",
        )
        image_rgba = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
        white = Image.new("RGBA", image_rgba.size, "#ffffff")
        white.alpha_composite(image_rgba)
        image = white.convert("RGB")

        image.thumbnail((640, 360), Image.Resampling.LANCZOS)
        thumbnail = Image.new("RGB", (640, 360), "#ffffff")
        left = (640 - image.width) // 2
        top = (360 - image.height) // 2
        thumbnail.paste(image, (left, top))
    except Exception:
        thumbnail = Image.new("RGB", (640, 360), "#ffffff")
        draw = ImageDraw.Draw(thumbnail)
        draw.rectangle((24, 28, 616, 332), outline="#b9c7ca", width=2)
        draw.text((40, 44), example.title, fill="#17202a", font=ImageFont.load_default())
        draw.text((40, 78), ", ".join(example.tags), fill="#1d6f7a", font=ImageFont.load_default())
    thumbnail.save(example.thumbnail_path, "WEBP", quality=82, method=6)


def _validate_source_assets(example: GallerySessionExample) -> None:
    missing = [
        str(path.relative_to(REPO_ROOT))
        for path in (example.session_path, example.output_svg_path)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing gallery source asset(s): {', '.join(missing)}")


def prepare_gallery_assets() -> list[dict[str, object]]:
    EXAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
    SESSION_ROOT.mkdir(parents=True, exist_ok=True)
    SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
    THUMBNAIL_ROOT.mkdir(parents=True, exist_ok=True)
    _remove_stale_assets()

    payload: list[dict[str, object]] = []
    for example in EXAMPLES:
        _validate_source_assets(example)
        session = _load_session(example)
        _render_thumbnail(example)

        source_figure = (
            str(example.source_svg_path.relative_to(REPO_ROOT))
            if example.source_svg_path.exists()
            else ""
        )
        payload.append(
            {
                "id": example.id,
                "title": example.title,
                "tags": list(example.tags),
                "description": example.description,
                "svg": example.gallery_svg_ref,
                "session": example.session_ref,
                "thumbnail": example.thumbnail_ref,
                "sourceSession": str(example.session_path.relative_to(REPO_ROOT)),
                "sourceOutput": str(example.output_svg_path.relative_to(REPO_ROOT)),
                "sourceFigure": source_figure,
                "sourceNote": example.source_note,
                "featureSources": _session_feature_sources(session),
                "fileSizeLabel": _format_size(example.gallery_svg_path.stat().st_size),
                "command": _session_command(session),
                "interactiveStep": example.interactive_step,
            }
        )

    (GALLERY_ROOT / "examples.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return payload


def main() -> int:
    payload = prepare_gallery_assets()
    print(f"Prepared {len(payload)} interactive gallery examples in {GALLERY_ROOT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
