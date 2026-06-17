#!/usr/bin/env python3
from __future__ import annotations

import base64
import io
import json
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import cairosvg
from PIL import Image, ImageDraw, ImageFont


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.web_support.feature_metadata import extract_features_from_genbank_payload  # noqa: E402


SOURCE_ROOT = REPO_ROOT / "examples"
TEST_INPUT_ROOT = REPO_ROOT / "tests" / "test_inputs"
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
EXAMPLE_ROOT = GALLERY_ROOT / "examples"
THUMBNAIL_ROOT = GALLERY_ROOT / "thumbnails"
STANDALONE_INTERACTIVITY_MODULE = WEB_ROOT / "js" / "services" / "standalone-interactivity.js"

INPUT_ALIASES = {
    "NC_012920.gb": TEST_INPUT_ROOT / "HmmtDNA.gbk",
}
PATH_TOKEN_SUFFIXES = (
    ".gb",
    ".gbk",
    ".gbff",
    ".gff",
    ".gff3",
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".ffn",
    ".faa",
    ".tsv",
    ".txt",
    ".out",
)


@dataclass(frozen=True)
class GalleryExample:
    id: str
    title: str
    tags: tuple[str, ...]
    description: str
    source_svg: str
    command: str
    feature_sources: tuple[str, ...]
    selected_features: tuple[str, ...] = ()
    interactive_step: str = (
        "Generated with the same Rich Interactive SVG enrichment used by the web app export."
    )
    source_note: str = "Base SVG is regenerated from the command during gallery asset preparation."
    thumbnail_fit: str = "contain"

    @property
    def svg_path(self) -> str:
        return f"./examples/{self.id}.svg"

    @property
    def thumbnail_path(self) -> str:
        return f"./thumbnails/{self.id}.webp"


EXAMPLES: tuple[GalleryExample, ...] = (
    GalleryExample(
        id="lambda-phage-linear",
        title="Escherichia phage Lambda linear map",
        tags=("Linear", "Lightweight"),
        description="A compact linear genome map for quick hover, click, pan, and zoom testing.",
        source_svg="NC_001416.svg",
        feature_sources=("NC_001416.gb",),
        command=(
            "gbdraw linear --gbk NC_001416.gb --show_labels all --separate_strands "
            "--legend left -d cds_white.tsv -t lambda_specific_table.tsv "
            "--block_stroke_width 2 --axis_stroke_width 5 --definition_font_size 24 "
            "-f svg -o NC_001416"
        ),
    ),
    GalleryExample(
        id="human-mtdna-compact",
        title="Human mitochondrial genome, compact circular map",
        tags=("Circular", "Compact", "Labels"),
        description="A small, label-rich circular genome for inspecting organellar labels and feature classes.",
        source_svg="NC_012920_middle_qualifier_priority_inner_axis5_def28_italic.svg",
        feature_sources=("NC_012920.gb",),
        command=(
            "gbdraw circular --gbk NC_012920.gb -f svg --track_type middle "
            '--species "<i>Homo sapiens</i>" --block_stroke_width 2 '
            "--axis_stroke_width 5 --labels both --qualifier_priority qualifier_priority.tsv "
            "-o NC_012920_middle_qualifier_priority_inner_axis5_def28_italic "
            "--definition_font_size 28"
        ),
    ),
    GalleryExample(
        id="hepatoplasmataceae-comparison",
        title="Hepatoplasmataceae five-genome comparison",
        tags=("Linear", "Comparison", "Featured", "Large"),
        description="A five-genome comparison spanning AP027078, AP027131, AP027133, AP027132, and NZ_CP006932.",
        source_svg="hepatoplasmataceae_default.svg",
        feature_sources=("AP027078.gb", "AP027131.gb", "AP027133.gb", "AP027132.gb", "NZ_CP006932.gb"),
        command=(
            "gbdraw linear --gbk AP027078.gb AP027131.gb AP027133.gb AP027132.gb NZ_CP006932.gb "
            "-b AP027078_AP027131.tblastx.out AP027131_AP027133.tblastx.out "
            "AP027133_AP027132.tblastx.out AP027132_NZ_CP006932.tblastx.out "
            "--align_center --separate_strands --block_stroke_width 1 "
            "--block_stroke_color gray --palette default -f svg -o hepatoplasmataceae_default"
        ),
    ),
    GalleryExample(
        id="majanivirus-comparison",
        title="Majanivirus ten-genome comparison",
        tags=("Linear", "Comparison", "Large"),
        description="A dense ten-genome comparison with many match blocks for stress-testing interactive inspection.",
        source_svg="majani.svg",
        feature_sources=(
            "MjeNMV.gb",
            "MelaMJNV.gb",
            "PemoMJNVA.gb",
            "PeseMJNV.gb",
            "PemoMJNVB.gb",
            "LvMJNV.gb",
            "TrcuMJNV.gb",
            "MellatMJNV.gb",
            "MeenMJNV.gb",
            "MejoMJNV.gb",
        ),
        command=(
            "gbdraw linear --gbk MjeNMV.gb MelaMJNV.gb PemoMJNVA.gb PeseMJNV.gb PemoMJNVB.gb "
            "LvMJNV.gb TrcuMJNV.gb MellatMJNV.gb MeenMJNV.gb MejoMJNV.gb "
            "-b MjeNMV.MelaMJNV.tblastx.out MelaMJNV.PemoMJNVA.tblastx.out "
            "PemoMJNVA.PeseMJNV.tblastx.out PeseMJNV.PemoMJNVB.tblastx.out "
            "PemoMJNVB.LvMJNV.tblastx.out LvMJNV.TrcuMJNV.tblastx.out "
            "TrcuMJNV.MellatMJNV.tblastx.out MellatMJNV.MeenMJNV.tblastx.out "
            "MeenMJNV.MejoMJNV.tblastx.out -t majani_custom_color_table.tsv "
            "-d modified_default_colors.tsv --block_stroke_width 1 "
            "--block_stroke_color gray --align_center --separate_strands -f svg -o majani"
        ),
    ),
    GalleryExample(
        id="tobacco-chloroplast",
        title="Nicotiana tabacum chloroplast genome",
        tags=("Circular", "Organellar", "Feature types"),
        description="A chloroplast example with clear feature classes, custom colors, and compact circular layout.",
        source_svg="NC_001879_color.svg",
        feature_sources=("NC_001879.gbk",),
        selected_features=("CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "rep_origin"),
        command=(
            "gbdraw circular --gbk NC_001879.gbk --separate_strands -f svg -o NC_001879_color "
            "-k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin -t chloroplast_specific_table.tsv "
            "--block_stroke_width 1 --block_stroke_color black --axis_stroke_width 3 "
            "--line_stroke_width 2 --suppress_gc --suppress_skew -p default "
            "--track_type tuckin --labels both --qualifier_priority qualifier_priority.tsv "
            "--outer_label_x_radius_offset 0.90 --outer_label_y_radius_offset 0.90 "
            "--inner_label_x_radius_offset 0.975 --inner_label_y_radius_offset 0.975 "
            '--species "<i>Nicotiana tabacum</i>" --definition_font_size 28 --legend upper_left'
        ),
    ),
)


@dataclass
class GalleryHarness:
    module_source: str
    _playwright: Any = field(default=None, init=False, repr=False)
    _browser: Any = field(default=None, init=False, repr=False)
    _page: Any = field(default=None, init=False, repr=False)

    def __enter__(self) -> "GalleryHarness":
        try:
            from playwright.sync_api import sync_playwright
        except ImportError as exc:
            raise RuntimeError(
                "Playwright is required to build gallery SVG assets. "
                "Install it with `python -m pip install playwright` and "
                "`python -m playwright install chromium`."
            ) from exc

        self._playwright = sync_playwright().start()
        try:
            self._browser = self._playwright.chromium.launch()
        except Exception as exc:
            self._playwright.stop()
            raise RuntimeError(
                "Could not launch Playwright Chromium. Run `python -m playwright install chromium`."
            ) from exc
        self._page = self._browser.new_page()
        self._page.set_content("<!doctype html><html><body></body></html>")
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self._browser is not None:
            self._browser.close()
        if self._playwright is not None:
            self._playwright.stop()

    def enrich_svg(self, svg_source: str, features: list[dict[str, object]]) -> str:
        module_url = "data:text/javascript;base64," + base64.b64encode(
            self.module_source.encode("utf-8")
        ).decode("ascii")
        return self._page.evaluate(
            """
            async ({ moduleUrl, svgSource, features }) => {
              const module = await import(moduleUrl);
              const parser = new DOMParser();
              const doc = parser.parseFromString(svgSource, 'image/svg+xml');
              const parserError = doc.querySelector('parsererror');
              if (parserError) {
                throw new Error(parserError.textContent || 'Could not parse SVG source.');
              }
              const svg = doc.documentElement;
              module.enrichSvgWithStandaloneInteractivity(svg, {
                features,
                popupMode: 'rich'
              });
              return new XMLSerializer().serializeToString(svg);
            }
            """,
            {
                "moduleUrl": module_url,
                "svgSource": svg_source,
                "features": features,
            },
        )


def _format_size(num_bytes: int) -> str:
    if num_bytes >= 1024 * 1024:
        return f"{num_bytes / (1024 * 1024):.1f} MB"
    return f"{num_bytes / 1024:.0f} KB"


def _resolve_input_path(filename: str) -> Path:
    alias = INPUT_ALIASES.get(filename)
    candidates = [
        alias,
        SOURCE_ROOT / filename,
        TEST_INPUT_ROOT / filename,
    ]
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate gallery input file: {filename}")


def _resolve_command_path_token(token: str) -> str:
    if not token.lower().endswith(PATH_TOKEN_SUFFIXES):
        return token
    try:
        return str(_resolve_input_path(token))
    except FileNotFoundError:
        return token


def _build_render_args(example: GalleryExample, output_prefix: Path) -> list[str]:
    tokens = shlex.split(example.command)
    if tokens and tokens[0] == "gbdraw":
        tokens = tokens[1:]
    if not tokens:
        raise ValueError(f"Gallery example {example.id} has an empty command.")

    args: list[str] = []
    output_arg_seen = False
    i = 0
    while i < len(tokens):
        token = tokens[i]
        args.append(token)
        if token in {"-o", "--output"}:
            output_arg_seen = True
            i += 1
            if i >= len(tokens):
                raise ValueError(f"Gallery example {example.id} has {token} without a value.")
            args.append(str(output_prefix))
        else:
            args[-1] = _resolve_command_path_token(token)
        i += 1

    if not output_arg_seen:
        args.extend(["-o", str(output_prefix)])
    return args


def _generate_base_svg(example: GalleryExample) -> str:
    with tempfile.TemporaryDirectory(prefix=f"gbdraw-gallery-{example.id}-") as temp_dir:
        output_prefix = Path(temp_dir) / example.id
        args = _build_render_args(example, output_prefix)
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", *args],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            message = "\n".join(part for part in (result.stdout, result.stderr) if part)
            raise RuntimeError(f"Could not regenerate gallery base SVG for {example.id}.\n{message}")
        svg_paths = sorted(output_prefix.parent.glob(f"{output_prefix.name}*.svg"))
        if not svg_paths:
            message = "\n".join(part for part in (result.stdout, result.stderr) if part)
            raise RuntimeError(f"No regenerated SVG found for {example.id}.\n{message}")
        return svg_paths[0].read_text(encoding="utf-8")


def _extract_feature_payload(example: GalleryExample) -> list[dict[str, object]]:
    features: list[dict[str, object]] = []
    record_offset = 0
    selected_features = list(example.selected_features) if example.selected_features else None
    for source in example.feature_sources:
        payload = extract_features_from_genbank_payload(
            _resolve_input_path(source),
            selected_features=selected_features,
        )
        source_features = payload.get("features", [])
        if not isinstance(source_features, list):
            continue
        for feature in source_features:
            if not isinstance(feature, dict):
                continue
            enriched = dict(feature)
            try:
                enriched["record_idx"] = int(enriched.get("record_idx") or 0) + record_offset
            except Exception:
                enriched["record_idx"] = record_offset
            enriched["id"] = f"f{len(features)}"
            features.append(enriched)
        record_ids = payload.get("record_ids", [])
        record_offset += len(record_ids) if isinstance(record_ids, list) else 1
    return features


def _render_thumbnail(example: GalleryExample, svg_source: str) -> None:
    output_path = THUMBNAIL_ROOT / f"{example.id}.webp"
    try:
        png_bytes = cairosvg.svg2png(
            bytestring=svg_source.encode("utf-8"),
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
    thumbnail.save(output_path, "WEBP", quality=82, method=6)


def prepare_gallery_assets() -> list[dict[str, object]]:
    EXAMPLE_ROOT.mkdir(parents=True, exist_ok=True)
    THUMBNAIL_ROOT.mkdir(parents=True, exist_ok=True)
    module_source = STANDALONE_INTERACTIVITY_MODULE.read_text(encoding="utf-8")
    payload: list[dict[str, object]] = []
    with GalleryHarness(module_source) as harness:
        for example in EXAMPLES:
            output_svg_path = EXAMPLE_ROOT / f"{example.id}.svg"
            source_svg = _generate_base_svg(example)
            features = _extract_feature_payload(example)
            output_svg_path.write_text(
                harness.enrich_svg(source_svg, features),
                encoding="utf-8",
            )
            _render_thumbnail(example, source_svg)
            payload.append(
                {
                    "id": example.id,
                    "title": example.title,
                    "tags": list(example.tags),
                    "description": example.description,
                    "svg": example.svg_path,
                    "thumbnail": example.thumbnail_path,
                    "sourceFigure": f"examples/{example.source_svg}",
                    "sourceNote": example.source_note,
                    "featureSources": list(example.feature_sources),
                    "fileSizeLabel": _format_size(output_svg_path.stat().st_size),
                    "command": example.command,
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
