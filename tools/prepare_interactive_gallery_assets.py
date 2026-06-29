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

from gbdraw.render.interactive_svg import InteractiveSvgContext, enrich_svg


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
    description: str = ""
    interactive_step: str = ""
    source_note: str = "Session JSON and generated SVG output are stored with the gallery assets."
    generate_output_from_session: bool = False

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
        title="<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018",
        tags=("Circular", "Multi-record"),
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_collinear",
        title="<i>Mollicutes</i> (<i>Candidatus</i> Hepatoplasmataceae) (collinear analysis)",
        tags=("Linear", "LOSAT", "Collinear"),
    ),
    GallerySessionExample(
        id="hepatoplasmataceae_orthogroup",
        title="<i>Mollicutes</i> (<i>Candidatus</i> Hepatoplasmataceae) (orthogroup matches)",
        tags=("Linear", "LOSAT", "Orthogroup"),
    ),
    GallerySessionExample(
        id="majanivirus_orthogroup",
        title="Large dsDNA viruses (<i>Nimaviridae</i>)",
        tags=("Linear", "LOSAT", "Orthogroup"),
    ),
    GallerySessionExample(
        id="BGC0000708-BGC0000713",
        title="Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.",
        tags=("Linear", "LOSAT", "Orthogroup"),
        generate_output_from_session=True,
    ),
    GallerySessionExample(
        id="HmmtDNA_ATskew",
        title="<i>Homo sapiens</i> mitochondrion (AT skew)",
        tags=("Circular", "Mitochondrion"),
        generate_output_from_session=True,
    ),
    GallerySessionExample(
        id="WSSV_genome_comparison",
        title="White spot syndrome virus genome comparison",
        tags=("Circular", "Conservation", "LOSAT"),
        generate_output_from_session=True,
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
    config = session.get("config") if isinstance(session.get("config"), dict) else {}
    ui = session.get("ui") if isinstance(session.get("ui"), dict) else {}
    cli_options = config.get("cliOptions") if isinstance(config.get("cliOptions"), dict) else {}
    mode = str(cli.get("mode") or cli_options.get("mode") or ui.get("mode") or "linear")
    raw_args = _session_raw_args(session)
    if raw_args:
        return shlex.join(["gbdraw", mode, *raw_args])
    return shlex.join(["gbdraw", mode, "--session", "session.gbdraw-session.json"])


def _session_feature_sources(session: dict[str, Any]) -> list[str]:
    cli = _session_cli_invocation(session)
    seen: set[str] = set()
    sources: list[str] = []

    def add_name(value: object) -> None:
        if not isinstance(value, str) or not value.lower().endswith(GENOME_SUFFIXES):
            return
        if value in seen:
            return
        seen.add(value)
        sources.append(value)

    bindings = cli.get("fileBindings")
    if isinstance(bindings, list):
        for binding in bindings:
            if not isinstance(binding, dict):
                continue
            add_name(binding.get("name"))

    files = session.get("files")
    if isinstance(files, dict):
        c_gb = files.get("c_gb")
        if isinstance(c_gb, dict):
            add_name(c_gb.get("name"))

        linear_seqs = files.get("linearSeqs")
        if isinstance(linear_seqs, list):
            for entry in linear_seqs:
                if not isinstance(entry, dict):
                    continue
                gb_entry = entry.get("gb")
                if isinstance(gb_entry, dict):
                    add_name(gb_entry.get("name"))
                add_name(entry.get("name"))
    return sources


def _session_result_svg(session: dict[str, Any], example: GallerySessionExample) -> str:
    results = session.get("results")
    if not isinstance(results, list):
        raise ValueError(f"{example.session_path} does not contain a results array.")
    for result in results:
        if not isinstance(result, dict):
            continue
        content = result.get("content")
        if isinstance(content, str) and "<svg" in content:
            return content
    raise ValueError(f"{example.session_path} does not contain a generated SVG result.")


def _session_interactive_context(session: dict[str, Any]) -> InteractiveSvgContext:
    feature_state = session.get("features") if isinstance(session.get("features"), dict) else {}
    editor_state = session.get("editorState") if isinstance(session.get("editorState"), dict) else {}
    orthogroup_state = (
        session.get("orthogroupState") if isinstance(session.get("orthogroupState"), dict) else {}
    )
    ui = session.get("ui") if isinstance(session.get("ui"), dict) else {}
    config = session.get("config") if isinstance(session.get("config"), dict) else {}
    legend = editor_state.get("legend") if isinstance(editor_state.get("legend"), dict) else {}

    features = feature_state.get("extractedFeatures")
    orthogroups = orthogroup_state.get("groups")
    legend_entries = legend.get("entries")
    current_colors = ui.get("appliedPaletteColors") or config.get("colors")
    return InteractiveSvgContext(
        features=features if isinstance(features, list) else (),
        orthogroups=orthogroups if isinstance(orthogroups, list) else (),
        legend_entries=legend_entries if isinstance(legend_entries, list) else (),
        current_colors=current_colors if isinstance(current_colors, dict) else {},
    )


def _write_source_svg(example: GallerySessionExample, session: dict[str, Any]) -> None:
    example.source_svg_path.write_text(_session_result_svg(session, example), encoding="utf-8")


def _write_gallery_svg_if_missing(example: GallerySessionExample, session: dict[str, Any]) -> None:
    if example.gallery_svg_path.exists() and not example.generate_output_from_session:
        return
    source = _session_result_svg(session, example)
    enriched = enrich_svg(source, context=_session_interactive_context(session))
    example.gallery_svg_path.write_text(enriched, encoding="utf-8")


def _remove_stale_assets() -> None:
    expected_svgs = {f"{example.id}.svg" for example in EXAMPLES}
    expected_thumbnails = {f"{example.id}.webp" for example in EXAMPLES}

    for path in EXAMPLE_ROOT.glob("*.svg"):
        if path.name not in expected_svgs:
            path.unlink()
    for path in SOURCE_ROOT.glob("*.svg"):
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
        for path in (example.session_path, example.output_svg_path, example.source_svg_path)
        if not path.exists()
    ]
    if missing:
        raise FileNotFoundError(f"Missing gallery source asset(s): {', '.join(missing)}")


def _gallery_entry_size(entry: dict[str, object]) -> int:
    svg_ref = str(entry["svg"])
    return (GALLERY_ROOT / svg_ref.removeprefix("./")).stat().st_size


def prepare_gallery_assets() -> list[dict[str, object]]:
    EXAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
    SESSION_ROOT.mkdir(parents=True, exist_ok=True)
    SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
    THUMBNAIL_ROOT.mkdir(parents=True, exist_ok=True)
    _remove_stale_assets()

    payload: list[dict[str, object]] = []
    for example in EXAMPLES:
        session = _load_session(example)
        _write_source_svg(example, session)
        _write_gallery_svg_if_missing(example, session)
        _validate_source_assets(example)
        _render_thumbnail(example)

        entry = {
            "id": example.id,
            "title": example.title,
            "tags": list(example.tags),
            "svg": example.gallery_svg_ref,
            "session": example.session_ref,
            "thumbnail": example.thumbnail_ref,
            "sourceSession": str(example.session_path.relative_to(REPO_ROOT)),
            "sourceOutput": str(example.output_svg_path.relative_to(REPO_ROOT)),
            "sourceFigure": str(example.source_svg_path.relative_to(REPO_ROOT)),
            "sourceNote": example.source_note,
            "featureSources": _session_feature_sources(session),
            "fileSizeLabel": _format_size(example.gallery_svg_path.stat().st_size),
            "command": _session_command(session),
        }
        if example.description:
            entry["description"] = example.description
        if example.interactive_step:
            entry["interactiveStep"] = example.interactive_step
        payload.append(entry)

    payload.sort(key=_gallery_entry_size)

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
