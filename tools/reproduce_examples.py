from __future__ import annotations

import argparse
import binascii
import hashlib
import inspect
import json
import math
import shlex
import struct
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from types import ModuleType
from typing import Any, Callable

from PIL import Image

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools import test_data  # noqa: E402
from tools.reproduce_examples_manifest import (  # noqa: E402
    CliRecipe,
    CompositePanel,
    CompositeRecipe,
    FigureSpec,
    SessionVariantRecipe,
    build_figure_specs,
    get_support_assets,
)


Image.MAX_IMAGE_PIXELS = None
RESAMPLING = getattr(Image, "Resampling", Image)

ALIASES = {
    "NC_000913.gbk": PROJECT_ROOT / "tests" / "test_inputs" / "MG1655.gbk",
    "Escherichia_coli.gbk": PROJECT_ROOT / "tests" / "test_inputs" / "MG1655.gbk",
    "O157_H7.gbk": PROJECT_ROOT / "tests" / "test_inputs" / "Sakai.gbk",
    "NC_012920.gb": PROJECT_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk",
}


# The PNG text chunk that names what a generated PNG was rasterized from. The pixels depend on the
# machine's FreeType and fonts, so `--check` compares this digest for PNGs, never their bytes.
PNG_SOURCE_KEY = "gbdraw-example-source"
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def source_digest(parts: list[bytes]) -> str:
    """sha256 over length-prefixed parts, so no two part lists share a digest."""
    digest = hashlib.sha256()
    for part in parts:
        digest.update(struct.pack(">Q", len(part)))
        digest.update(part)
    return digest.hexdigest()


def png_with_source_digest(png: bytes, digest: str) -> bytes:
    """``png`` with a ``PNG_SOURCE_KEY`` tEXt chunk right after IHDR, before the image data, so it
    is read without decoding the image."""
    if png[:8] != _PNG_SIGNATURE or png[12:16] != b"IHDR":
        raise ValueError("Not a PNG.")
    data = PNG_SOURCE_KEY.encode("latin-1") + b"\0" + digest.encode("ascii")
    chunk = struct.pack(">I", len(data)) + b"tEXt" + data + struct.pack(">I", binascii.crc32(b"tEXt" + data))
    end_of_ihdr = 8 + 12 + struct.unpack(">I", png[8:12])[0]
    return png[:end_of_ihdr] + chunk + png[end_of_ihdr:]


def png_source_digest(png: bytes) -> str | None:
    """The ``PNG_SOURCE_KEY`` digest a PNG carries before its image data, or None."""
    if png[:8] != _PNG_SIGNATURE:
        return None
    offset = 8
    prefix = PNG_SOURCE_KEY.encode("latin-1") + b"\0"
    while offset + 8 <= len(png):
        length, kind = struct.unpack(">I4s", png[offset : offset + 8])
        if kind == b"IDAT":
            return None
        data = png[offset + 8 : offset + 8 + length]
        if kind == b"tEXt" and data.startswith(prefix):
            return data[len(prefix) :].decode("ascii")
        offset += 12 + length
    return None


def raster_input_digest(options: dict[str, Any]) -> bytes:
    """What one CairoSVG PNG conversion draws from, as gbdraw calls ``svg2png``: the SVG bytes and the
    options that size the image (scale, output size). ``url`` only locates relative references."""
    sizing = sorted((key, value) for key, value in options.items() if key not in {"bytestring", "url"})
    return source_digest([options["bytestring"], repr(sizing).encode("utf-8")]).encode("ascii")


def recording_svg2png(svg2png: Callable[..., Any], record: Callable[[bytes], None]) -> Callable[..., Any]:
    """``svg2png`` that first hands ``record`` the ``raster_input_digest`` of each conversion: the one
    record of a PNG's raster input, whatever gbdraw code prepared the bytes."""

    def convert(**options: Any) -> Any:
        record(raster_input_digest(options))
        return svg2png(**options)

    return convert


def run_recording_cli(argv: list[str]) -> None:
    """``gbdraw <argv[1:]>`` with the ``raster_input_digest`` of each PNG it writes appended, one per
    line, to the file ``argv[0]``."""
    from gbdraw.cli import main
    from gbdraw.render.export import get_cairosvg

    record = Path(argv[0])

    def append(digest: bytes) -> None:
        with record.open("ab") as handle:
            handle.write(digest + b"\n")

    # The module gbdraw's export converts with: the wrapper sees each call, and converts nothing itself.
    cairosvg_module = get_cairosvg()
    cairosvg_module.svg2png = recording_svg2png(cairosvg_module.svg2png, append)
    sys.argv = ["gbdraw", *argv[1:]]
    main()


# ``python -c RECORDING_CLI <record file> <subcommand> ...`` from the project root: a CLI PNG recipe.
RECORDING_CLI = "import sys; from tools.reproduce_examples import run_recording_cli; run_recording_cli(sys.argv[1:])"


def tracked_input_path(project_root: Path, filename: str) -> Path | None:
    """An input file the repository (or the test-data Release) places under examples/ or
    tests/test_inputs/."""
    for directory in (project_root / "examples", project_root / "tests" / "test_inputs"):
        candidate = directory / filename
        if candidate.exists():
            return candidate
    return None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce README/docs/palette example figures for gbdraw.",
    )
    parser.add_argument(
        "--group",
        action="append",
        choices=("docs", "palettes", "composites", "readme", "all"),
        default=None,
        help="Figure group to render. Repeatable. Defaults to all.",
    )
    parser.add_argument(
        "--figure",
        action="append",
        default=None,
        help="Specific figure ID to render. Repeatable.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Root directory for reproduced outputs (default: _reproduced, or a temporary "
            "directory with --check). '.' writes the tracked files in place."
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Render the figures and compare each with the tracked file (an SVG byte for byte, a "
            "PNG by the digest of what it was rasterized from); list every figure that differs "
            "and the command that refreshes it, and exit 1."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Figures rendered at once, each in its own process (default: %(default)s).",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List matching figure IDs and readiness without rendering.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip rendering when the destination output already exists.",
    )
    parser.add_argument(
        "--missing-report",
        type=Path,
        default=None,
        help="Path to write the missing-input JSON report.",
    )
    return parser.parse_args(argv)


class Reproducer:
    def __init__(
        self,
        project_root: Path,
        output_root: Path,
        figures: dict[str, FigureSpec],
        skip_existing: bool = False,
    ) -> None:
        self.project_root = project_root
        self.output_root = output_root
        self.figures = figures
        self.skip_existing = skip_existing
        self.support_assets = get_support_assets()
        self.generated: list[str] = []
        self.skipped_missing_inputs: list[str] = []
        self.failed: list[str] = []
        self.aliases_used: dict[tuple[str, str], set[str]] = defaultdict(set)
        self.missing_inputs: dict[str, dict[str, Any]] = {}
        self._render_state: dict[str, str] = {}
        self._inspect_state: dict[str, bool] = {}
        self._active_render: set[str] = set()
        self._active_inspect: set[str] = set()
        self._tmpdir = tempfile.TemporaryDirectory(prefix="reproduce_examples_")
        self.temp_root = Path(self._tmpdir.name)
        # Intermediate inputs stay out of the output root, which may be the repository itself.
        self.support_dir = self.temp_root / "_support_assets"

    def close(self) -> None:
        self._tmpdir.cleanup()

    def mark_rendered(self, figure_ids: set[str]) -> None:
        """Treat figures another process already rendered into the output root as done."""
        for figure_id in figure_ids:
            self._render_state[figure_id] = "existing"

    def output_path_for(self, figure_id: str) -> Path:
        return self.output_root / self.figures[figure_id].output_path

    def report_payload(self) -> dict[str, Any]:
        alias_entries = [
            {
                "requested": requested,
                "resolved_to": resolved,
                "figures": sorted(figures),
            }
            for (requested, resolved), figures in sorted(self.aliases_used.items())
        ]
        missing_entries = []
        for filename, data in sorted(self.missing_inputs.items()):
            missing_entries.append(
                {
                    "filename": filename,
                    "figures": sorted(data["figures"]),
                }
            )
        return {
            "generated": self.generated,
            "skipped_missing_inputs": self.skipped_missing_inputs,
            "failed": self.failed,
            "aliases_used": alias_entries,
            "missing_inputs": missing_entries,
        }

    def merge_report(self, payload: dict[str, Any]) -> None:
        """Fold the report of a Reproducer that rendered other figures into this one."""
        for key, target in (
            ("generated", self.generated),
            ("skipped_missing_inputs", self.skipped_missing_inputs),
            ("failed", self.failed),
        ):
            target.extend(figure_id for figure_id in payload[key] if figure_id not in target)
        for alias in payload["aliases_used"]:
            self.aliases_used[(alias["requested"], alias["resolved_to"])].update(alias["figures"])
        for missing in payload["missing_inputs"]:
            for figure_id in missing["figures"]:
                self.add_missing(missing["filename"], figure_id)

    def write_report(self, report_path: Path) -> None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(self.report_payload(), indent=2, sort_keys=False) + "\n")

    def add_missing(self, filename: str, figure_id: str) -> None:
        self.missing_inputs.setdefault(filename, {"figures": set()})["figures"].add(figure_id)

    def add_alias_use(self, requested: str, resolved: Path, figure_id: str) -> None:
        resolved_str = str(resolved.relative_to(self.project_root))
        self.aliases_used[(requested, resolved_str)].add(figure_id)

    def _existing_input_path(self, filename: str) -> Path | None:
        tracked = tracked_input_path(self.project_root, filename)
        if tracked is not None:
            return tracked
        candidate = self.support_dir / filename
        return candidate if candidate.exists() else None

    def _alias_input_path(self, filename: str, figure_id: str) -> Path | None:
        alias_target = ALIASES.get(filename)
        if alias_target and alias_target.exists():
            self.add_alias_use(filename, alias_target, figure_id)
            return alias_target
        return None

    def _materialize_support_asset(self, filename: str, dry_run: bool) -> Path | None:
        asset = self.support_assets.get(filename)
        if not asset:
            return None
        target_path = self.support_dir / asset.filename
        if dry_run:
            return target_path
        target_path.parent.mkdir(parents=True, exist_ok=True)
        content = asset.content.replace("{PROJECT_ROOT}", str(self.project_root))
        if not target_path.exists() or target_path.read_text(encoding="utf-8") != content:
            target_path.write_text(content, encoding="utf-8")
        return target_path

    def resolve_input(self, filename: str, figure_id: str, *, dry_run: bool) -> Path | None:
        asset_path = self._materialize_support_asset(filename, dry_run)
        if asset_path:
            return asset_path

        existing = self._existing_input_path(filename)
        if existing:
            return existing

        alias = self._alias_input_path(filename, figure_id)
        if alias:
            return alias

        self.add_missing(filename, figure_id)
        return None

    def _resolve_cli_recipe(
        self,
        recipe: CliRecipe,
        figure_id: str,
        *,
        dry_run: bool,
    ) -> dict[str, Any] | None:
        resolved: dict[str, Any] = {
            "gbk_files": [],
            "gff_files": [],
            "fasta_files": [],
            "blast_files": [],
            "file_args": [],
            "extra_args": list(recipe.extra_args),
        }
        categories = (
            ("gbk_files", recipe.gbk_files),
            ("gff_files", recipe.gff_files),
            ("fasta_files", recipe.fasta_files),
            ("blast_files", recipe.blast_files),
        )
        for key, filenames in categories:
            for filename in filenames:
                path = self.resolve_input(filename, figure_id, dry_run=dry_run)
                if path is None:
                    return None
                resolved[key].append(path)
        for file_arg in recipe.file_args:
            resolved_paths: list[Path] = []
            for filename in file_arg.filenames:
                path = self.resolve_input(filename, figure_id, dry_run=dry_run)
                if path is None:
                    return None
                resolved_paths.append(path)
            resolved["file_args"].append((file_arg.flag, resolved_paths))
        return resolved

    def _build_command(
        self,
        recipe: CliRecipe,
        resolved: dict[str, Any],
        output_prefix: Path,
        fmt: str,
    ) -> list[str]:
        if fmt == "png":
            launcher = ["-c", RECORDING_CLI, str(output_prefix.with_suffix(".raster"))]
        else:
            launcher = ["-m", "gbdraw.cli"]
        cmd = [sys.executable, *launcher, recipe.subcommand]
        if resolved["gbk_files"]:
            cmd.extend(["--gbk", *[str(path) for path in resolved["gbk_files"]]])
        if resolved["gff_files"]:
            cmd.extend(["--gff", *[str(path) for path in resolved["gff_files"]]])
        if resolved["fasta_files"]:
            cmd.extend(["--fasta", *[str(path) for path in resolved["fasta_files"]]])
        if resolved["blast_files"]:
            cmd.extend(["-b", *[str(path) for path in resolved["blast_files"]]])
        for flag, paths in resolved["file_args"]:
            cmd.extend([flag, *[str(path) for path in paths]])
        # A rerun into the same root, or a refresh in place, replaces the earlier output.
        cmd.extend(["-o", str(output_prefix), "-f", fmt, "--overwrite"])
        cmd.extend(resolved["extra_args"])
        return cmd

    def _run_cli_recipe(
        self,
        recipe: CliRecipe,
        figure_id: str,
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        resolved = self._resolve_cli_recipe(recipe, figure_id, dry_run=dry_run)
        if resolved is None:
            return False
        if dry_run:
            return True

        target_path.parent.mkdir(parents=True, exist_ok=True)
        fmt = target_path.suffix.lstrip(".")
        # A PNG renders in the temporary root, beside the CLI's SVG and the record of what CairoSVG
        # drew it from, and reaches the target with that record's digest.
        output_prefix = (
            self.temp_root / "cli_png" / figure_id / target_path.stem
            if fmt == "png"
            else target_path.with_suffix("")
        )
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        output_prefix.with_suffix(".raster").unlink(missing_ok=True)
        cmd = self._build_command(recipe, resolved, output_prefix, fmt)
        result = subprocess.run(
            cmd,
            cwd=str(self.project_root),
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or result.stdout.strip() or f"gbdraw {recipe.subcommand} failed")
        if fmt == "png":
            digest = source_digest(
                [
                    b"cli-png",
                    repr(recipe).encode("utf-8"),
                    *output_prefix.with_suffix(".raster").read_bytes().splitlines(),
                ]
            )
            target_path.write_bytes(
                png_with_source_digest(output_prefix.with_suffix(".png").read_bytes(), digest)
            )
        return True

    def _run_session_variant_recipe(
        self,
        recipe: SessionVariantRecipe,
        figure_id: str,
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        session_path = self.project_root / recipe.source_session_path
        if not session_path.is_file():
            self.add_missing(recipe.source_session_path, figure_id)
            return False
        if target_path.suffix.lower() != ".svg":
            raise RuntimeError("Session variant recipes require an SVG target")
        if dry_run:
            return True

        session = json.loads(session_path.read_text(encoding="utf-8"))
        # A Session 51 holds the request in its drawing; the variant edits it in place.
        active_id = session.get("ui", {}).get("activeDrawingId")
        drawings = [d for d in session.get("drawings", []) if d.get("id") == active_id]
        request = drawings[0].get("renderRequest") if len(drawings) == 1 else None
        resources = session.get("resources")
        if not isinstance(request, dict) or not isinstance(resources, dict):
            raise RuntimeError(
                f"Gallery session lacks a canonical request or resources: {session_path}"
            )
        diagram_options = request.get("diagramOptions")
        if not isinstance(diagram_options, dict):
            raise RuntimeError(
                f"Gallery session lacks diagramOptions: {session_path}"
            )
        feature_shapes = diagram_options.setdefault("featureShapes", {})
        config_overrides = diagram_options.setdefault("configOverrides", {})
        if not isinstance(feature_shapes, dict) or not isinstance(
            config_overrides, dict
        ):
            raise RuntimeError(
                f"Gallery session has invalid diagram options: {session_path}"
            )
        feature_shapes.update(recipe.feature_shapes)
        config_overrides.update(recipe.config_overrides)

        mode = request.get("mode")
        if mode not in {"circular", "linear"}:
            raise RuntimeError(f"Gallery session has an invalid mode: {session_path}")
        variant_session_path = (
            self.temp_root / f"{figure_id}.gbdraw-session.json"
        )
        variant_session_path.write_text(
            json.dumps(session, separators=(",", ":")),
            encoding="utf-8",
        )
        target_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable,
            "-m",
            "gbdraw.cli",
            mode,
            "--session",
            str(variant_session_path),
            "-o",
            str(target_path.with_suffix("")),
            "-f",
            "svg",
            "--overwrite",
        ]
        result = subprocess.run(
            cmd,
            cwd=str(self.project_root),
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            raise RuntimeError(
                result.stderr.strip()
                or result.stdout.strip()
                or f"Session variant render failed: {session_path}"
            )
        if not target_path.is_file():
            raise RuntimeError(f"Session variant did not render an SVG: {session_path}")
        return True

    def _svg_size(self, svg_path: Path) -> tuple[int, int]:
        import re

        text = svg_path.read_text(encoding="utf-8", errors="ignore")
        view_box_match = re.search(r'viewBox="([\d.\-]+) ([\d.\-]+) ([\d.\-]+) ([\d.\-]+)"', text)
        if view_box_match:
            width = int(round(float(view_box_match.group(3))))
            height = int(round(float(view_box_match.group(4))))
            return max(width, 1), max(height, 1)
        width_match = re.search(r'width="([\d.]+)(?:px)?"', text)
        height_match = re.search(r'height="([\d.]+)(?:px)?"', text)
        if width_match and height_match:
            width = int(round(float(width_match.group(1))))
            height = int(round(float(height_match.group(1))))
            return max(width, 1), max(height, 1)
        return (1000, 1000)

    def _load_image(
        self, source_path: Path, target_size: tuple[int, int], sources: list[bytes]
    ) -> Image.Image:
        """A panel image fitted to ``target_size``; appends what it was rasterized from to ``sources``."""
        sources.append(repr(target_size).encode("ascii"))
        if source_path.suffix.lower() == ".svg":
            from gbdraw.render.export import convert_svg_with_cairosvg, get_cairosvg

            try:
                cairosvg_module = get_cairosvg()
            except ImportError as exc:  # pragma: no cover - environment dependent
                raise RuntimeError("CairoSVG is required to build composite PNG outputs") from exc
            src_width, src_height = self._svg_size(source_path)
            scale = min(target_size[0] / src_width, target_size[1] / src_height)
            render_width = max(1, int(round(src_width * scale)))
            render_height = max(1, int(round(src_height * scale)))
            recording = ModuleType(cairosvg_module.__name__)
            recording.svg2png = recording_svg2png(cairosvg_module.svg2png, sources.append)
            png_bytes = convert_svg_with_cairosvg(
                source_path.read_bytes(),
                "png",
                cairosvg_module=recording,
                url=str(source_path),
                output_width=render_width,
                output_height=render_height,
            )
            image = Image.open(BytesIO(png_bytes)).convert("RGBA")
        else:
            raster = source_path.read_bytes()
            sources.append((png_source_digest(raster) or hashlib.sha256(raster).hexdigest()).encode("ascii"))
            image = Image.open(BytesIO(raster)).convert("RGBA")
        fitted = image.copy()
        fitted.thumbnail(target_size, RESAMPLING.LANCZOS)
        canvas = Image.new("RGBA", target_size, (255, 255, 255, 0))
        left = (target_size[0] - fitted.width) // 2
        top = (target_size[1] - fitted.height) // 2
        canvas.paste(fitted, (left, top), fitted)
        return canvas

    def _render_panel_source(
        self,
        panel: CompositePanel,
        figure_id: str,
        panel_index: int,
        target_size: tuple[int, int],
        sources: list[bytes],
        *,
        dry_run: bool,
    ) -> Image.Image | None:
        if panel.figure_id:
            if dry_run:
                return None if not self.inspect_figure(panel.figure_id) else Image.new("RGBA", target_size, (255, 255, 255, 0))
            if not self.render_figure(panel.figure_id):
                return None
            return self._load_image(self.output_path_for(panel.figure_id), target_size, sources)
        if panel.recipe is None:
            raise RuntimeError(f"Composite panel {panel_index} in {figure_id} has no source")
        if dry_run:
            return None if not self._run_cli_recipe(panel.recipe, figure_id, self.temp_root / "dry_run.svg", dry_run=True) else Image.new("RGBA", target_size, (255, 255, 255, 0))
        panel_dir = self.temp_root / figure_id
        panel_dir.mkdir(parents=True, exist_ok=True)
        panel_path = panel_dir / f"panel_{panel_index}.svg"
        if not self._run_cli_recipe(panel.recipe, figure_id, panel_path, dry_run=False):
            return None
        return self._load_image(panel_path, target_size, sources)

    def _render_grid_composite(
        self,
        figure_id: str,
        recipe: CompositeRecipe,
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        if recipe.columns is None or recipe.tile_size is None:
            raise RuntimeError(f"Composite recipe for {figure_id} is missing grid settings")
        if dry_run:
            for index, panel in enumerate(recipe.panels):
                image = self._render_panel_source(panel, figure_id, index, recipe.tile_size, [], dry_run=True)
                if image is None:
                    return False
            return True

        panel_images: list[Image.Image] = []
        sources: list[bytes] = []
        for index, panel in enumerate(recipe.panels):
            image = self._render_panel_source(panel, figure_id, index, recipe.tile_size, sources, dry_run=False)
            if image is None:
                return False
            panel_images.append(image)
        rows = int(math.ceil(len(panel_images) / recipe.columns))
        if recipe.canvas_size is None:
            width = recipe.padding * 2 + recipe.columns * recipe.tile_size[0] + (recipe.columns - 1) * recipe.gap
            height = recipe.padding * 2 + rows * recipe.tile_size[1] + (rows - 1) * recipe.gap
            canvas_size = (width, height)
        else:
            canvas_size = recipe.canvas_size
        background = Image.new("RGBA", canvas_size, recipe.background)
        for index, image in enumerate(panel_images):
            row = index // recipe.columns
            col = index % recipe.columns
            left = recipe.padding + col * (recipe.tile_size[0] + recipe.gap)
            top = recipe.padding + row * (recipe.tile_size[1] + recipe.gap)
            background.paste(image, (left, top), image)
        self._save_composite(recipe, background, target_path, sources)
        return True

    def _render_social_preview(
        self,
        figure_id: str,
        recipe: CompositeRecipe,
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        if recipe.canvas_size is None:
            raise RuntimeError("social_preview recipes require canvas_size")
        if dry_run:
            return all(
                panel.figure_id is not None and self.inspect_figure(panel.figure_id)
                for panel in recipe.panels
            )
        canvas = Image.new("RGBA", recipe.canvas_size, recipe.background)
        sources: list[bytes] = []
        for index, panel in enumerate(recipe.panels):
            if panel.figure_id is None or panel.box is None:
                raise RuntimeError(f"social_preview panel {index} for {figure_id} must reference a figure and a box")
            if not self.render_figure(panel.figure_id):
                return False
            box_width = panel.box[2]
            box_height = panel.box[3]
            image = self._load_image(self.output_path_for(panel.figure_id), (box_width, box_height), sources)
            canvas.paste(image, (panel.box[0], panel.box[1]), image)
        self._save_composite(recipe, canvas, target_path, sources)
        return True

    def _save_composite(
        self, recipe: CompositeRecipe, image: Image.Image, target_path: Path, sources: list[bytes]
    ) -> None:
        """Write a composite PNG with the digest of its recipe, the code that builds composites,
        and what each panel was rasterized from."""
        code = "".join(inspect.getsource(function) for function in _composite_code())
        digest = source_digest(
            [b"composite", repr(recipe).encode("utf-8"), code.encode("utf-8"), *sources]
        )
        buffer = BytesIO()
        image.convert("RGB").save(buffer, format="PNG")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        target_path.write_bytes(png_with_source_digest(buffer.getvalue(), digest))

    def _render_composite(
        self,
        figure_id: str,
        recipe: CompositeRecipe,
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        if recipe.kind in {"grid", "contact_sheet"}:
            return self._render_grid_composite(figure_id, recipe, target_path, dry_run=dry_run)
        if recipe.kind == "social_preview":
            return self._render_social_preview(figure_id, recipe, target_path, dry_run=dry_run)
        raise RuntimeError(f"Unsupported composite recipe kind: {recipe.kind}")

    def inspect_figure(self, figure_id: str) -> bool:
        if figure_id in self._inspect_state:
            return self._inspect_state[figure_id]
        if figure_id in self._active_inspect:
            raise RuntimeError(f"Cyclic figure dependency detected while inspecting {figure_id}")
        self._active_inspect.add(figure_id)
        spec = self.figures[figure_id]
        target_path = self.output_path_for(figure_id)
        try:
            if self.skip_existing and target_path.exists():
                ready = True
            elif isinstance(spec.recipe, CliRecipe):
                ready = self._run_cli_recipe(spec.recipe, figure_id, target_path, dry_run=True)
            elif isinstance(spec.recipe, SessionVariantRecipe):
                ready = self._run_session_variant_recipe(
                    spec.recipe,
                    figure_id,
                    target_path,
                    dry_run=True,
                )
            else:
                ready = self._render_composite(figure_id, spec.recipe, target_path, dry_run=True)
        finally:
            self._active_inspect.remove(figure_id)
        self._inspect_state[figure_id] = ready
        return ready

    def render_figure(self, figure_id: str) -> bool:
        current_state = self._render_state.get(figure_id)
        if current_state in {"generated", "existing"}:
            return True
        if current_state in {"skipped", "failed"}:
            return False
        if figure_id in self._active_render:
            raise RuntimeError(f"Cyclic figure dependency detected while rendering {figure_id}")

        self._active_render.add(figure_id)
        spec = self.figures[figure_id]
        target_path = self.output_path_for(figure_id)
        try:
            if self.skip_existing and target_path.exists():
                self._render_state[figure_id] = "existing"
                return True
            if isinstance(spec.recipe, CliRecipe):
                success = self._run_cli_recipe(spec.recipe, figure_id, target_path, dry_run=False)
            elif isinstance(spec.recipe, SessionVariantRecipe):
                success = self._run_session_variant_recipe(
                    spec.recipe,
                    figure_id,
                    target_path,
                    dry_run=False,
                )
            else:
                success = self._render_composite(figure_id, spec.recipe, target_path, dry_run=False)
            if not success:
                if figure_id not in self.skipped_missing_inputs:
                    self.skipped_missing_inputs.append(figure_id)
                self._render_state[figure_id] = "skipped"
                return False
            self.generated.append(figure_id)
            self._render_state[figure_id] = "generated"
            return True
        except Exception:
            if figure_id not in self.failed:
                self.failed.append(figure_id)
            self._render_state[figure_id] = "failed"
            raise
        finally:
            self._active_render.remove(figure_id)


def select_figures(figures: dict[str, FigureSpec], groups: list[str] | None, figure_ids: list[str] | None) -> list[str]:
    if groups is None or "all" in groups:
        selected = list(figures.keys())
    else:
        selected = [
            figure_id
            for figure_id, spec in figures.items()
            if any(group in spec.groups for group in groups)
        ]
    if figure_ids:
        missing_ids = sorted(set(figure_ids) - set(figures))
        if missing_ids:
            raise SystemExit(f"Unknown figure ID(s): {', '.join(missing_ids)}")
        selected = [figure_id for figure_id in selected if figure_id in set(figure_ids)]
    return selected


def print_list_status(reproducer: Reproducer, figure_ids: list[str]) -> None:
    for figure_id in figure_ids:
        output_path = reproducer.output_path_for(figure_id)
        if output_path.exists():
            status = "existing"
        else:
            status = "ready" if reproducer.inspect_figure(figure_id) else "missing"
        print(f"{figure_id}\t{status}\t{reproducer.figures[figure_id].output_path}")


def _composite_code() -> tuple[Any, ...]:
    """The functions whose code decides a composite's pixels beyond its panels' SVGs."""
    return (
        Reproducer._svg_size,
        Reproducer._load_image,
        Reproducer._render_panel_source,
        Reproducer._render_grid_composite,
        Reproducer._render_social_preview,
        Reproducer._save_composite,
    )


@dataclass(frozen=True)
class FigureRun:
    """What one process reports after rendering one figure (and the figures it depends on)."""

    figure_id: str
    report: dict[str, Any]
    error: str | None
    seconds: float


def _render_one(
    figures: dict[str, FigureSpec],
    figure_id: str,
    output_root: Path,
    rendered: frozenset[str],
    skip_existing: bool,
) -> FigureRun:
    reproducer = Reproducer(
        project_root=PROJECT_ROOT,
        output_root=output_root,
        figures=figures,
        skip_existing=skip_existing,
    )
    started = time.monotonic()
    error = None
    try:
        reproducer.mark_rendered(set(rendered))
        reproducer.render_figure(figure_id)
    except Exception as exc:
        error = str(exc) or type(exc).__name__
    finally:
        reproducer.close()
    return FigureRun(figure_id, reproducer.report_payload(), error, time.monotonic() - started)


class _InlineExecutor:
    """Runs each submitted call at once, so one process takes the same path as a pool."""

    def __enter__(self) -> _InlineExecutor:
        return self

    def __exit__(self, *exc_info: object) -> None:
        return None

    def submit(self, fn: Any, /, *args: Any) -> Future[FigureRun]:
        future: Future[FigureRun] = Future()
        future.set_result(fn(*args))
        return future


def _input_bytes(spec: FigureSpec) -> int:
    """The bytes of the inputs a figure reads, once per panel: the order the pool starts them in,
    so the figures that take longest do not start last."""
    size = 0
    for filename in spec.required_inputs:
        path = tracked_input_path(PROJECT_ROOT, filename) or ALIASES.get(filename)
        if path is not None and path.is_file():
            size += path.stat().st_size
    return size * (len(spec.recipe.panels) if isinstance(spec.recipe, CompositeRecipe) else 1)


def _panel_figures(spec: FigureSpec) -> set[str]:
    if not isinstance(spec.recipe, CompositeRecipe):
        return set()
    return {panel.figure_id for panel in spec.recipe.panels if panel.figure_id}


def render_figures(
    figures: dict[str, FigureSpec],
    figure_ids: list[str],
    output_root: Path,
    *,
    jobs: int,
    skip_existing: bool = False,
) -> list[FigureRun]:
    """Render each figure in its own Reproducer, ``jobs`` at a time; a composite starts after
    the selected figures it is built from, and reuses their outputs."""
    selected = set(figure_ids)
    waiting = sorted(figure_ids, key=lambda figure_id: -_input_bytes(figures[figure_id]))
    running: dict[Future[FigureRun], str] = {}
    finished: set[str] = set()
    rendered: set[str] = set()
    runs: list[FigureRun] = []
    executor = ProcessPoolExecutor(max_workers=jobs) if jobs > 1 else _InlineExecutor()
    with executor:
        while waiting or running:
            for figure_id in [
                figure_id
                for figure_id in waiting
                if (_panel_figures(figures[figure_id]) & selected) <= finished
            ]:
                waiting.remove(figure_id)
                future = executor.submit(
                    _render_one, figures, figure_id, output_root, frozenset(rendered), skip_existing
                )
                running[future] = figure_id
            if not running:
                raise RuntimeError(f"Cyclic figure dependency among: {', '.join(waiting)}")
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for future in done:
                run = future.result()
                del running[future]
                finished.add(run.figure_id)
                rendered.update(run.report["generated"])
                runs.append(run)
                status = "failed" if run.error else (
                    "generated" if run.figure_id in run.report["generated"] else "not rendered"
                )
                print(f"{run.figure_id}: {status} in {run.seconds:.1f} s", flush=True)
    return runs


def refresh_command(figure_ids: list[str]) -> str:
    arguments = ["python", "tools/reproduce_examples.py", "--output-root", "."]
    for figure_id in figure_ids:
        arguments.extend(["--figure", figure_id])
    return shlex.join(arguments)


def _same_output(tracked: bytes, generated: bytes) -> bool:
    generated_digest = png_source_digest(generated)
    if generated_digest is not None:
        return png_source_digest(tracked) == generated_digest
    return tracked == generated


def check_outputs(
    figures: dict[str, FigureSpec],
    figure_ids: list[str],
    output_root: Path,
    reproducer: Reproducer,
    errors: dict[str, str],
    *,
    tracked_root: Path = PROJECT_ROOT,
) -> bool:
    """Compare every rendered figure with its tracked file and print the result: an SVG byte for
    byte, a PNG by the source digest it carries (a tracked PNG without one is stale)."""
    same: list[str] = []
    stale: list[str] = []
    problems: list[str] = []
    for figure_id in figure_ids:
        spec = figures[figure_id]
        if figure_id in errors:
            problems.append(f"{figure_id}: rendering failed: {errors[figure_id]}")
        elif figure_id not in reproducer.generated:
            missing = sorted(
                filename
                for filename, entry in reproducer.missing_inputs.items()
                if figure_id in entry["figures"]
            )
            problems.append(f"{figure_id}: not rendered; missing inputs: {', '.join(missing) or 'none reported'}")
        else:
            tracked = tracked_root / spec.output_path
            generated = output_root / spec.output_path
            if tracked.is_file() and _same_output(tracked.read_bytes(), generated.read_bytes()):
                same.append(figure_id)
            else:
                stale.append(figure_id)
    print(
        f"\nExamples check: {len(same)} match, {len(stale)} differ, "
        f"{len(problems)} failed."
    )
    if stale:
        print("\nThese tracked files differ from what the generator renders now:")
        for figure_id in stale:
            tracked = tracked_root / figures[figure_id].output_path
            state = "" if tracked.is_file() else " (not tracked yet)"
            print(f"  {figures[figure_id].output_path}{state}")
        print(
            "\nReview the rendered figures, then refresh the tracked files with:\n  "
            + refresh_command(stale)
        )
    if problems:
        print("\nThese figures could not be checked:")
        for line in problems:
            print(f"  {line}")
    return not stale and not problems


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.check and (args.list or args.skip_existing):
        raise SystemExit("--check renders every selected figure; it takes neither --list nor --skip-existing.")
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1.")
    # Large inputs live in the test-data Release; a missing one is reported per figure below.
    test_data.fetch()
    figures = build_figure_specs()
    selected_ids = select_figures(figures, args.group, args.figure)
    with tempfile.TemporaryDirectory(prefix="reproduce_examples_check_") as check_root:
        if args.output_root is None:
            output_root = Path(check_root) if args.check else PROJECT_ROOT / "_reproduced"
        elif args.output_root.is_absolute():
            output_root = args.output_root
        else:
            output_root = PROJECT_ROOT / args.output_root
        missing_report = args.missing_report
        if missing_report is None:
            missing_report = output_root / "missing_inputs.json"
        elif not missing_report.is_absolute():
            missing_report = PROJECT_ROOT / missing_report

        reproducer = Reproducer(
            project_root=PROJECT_ROOT,
            output_root=output_root,
            figures=figures,
            skip_existing=args.skip_existing,
        )
        try:
            if args.list:
                print_list_status(reproducer, selected_ids)
                ok = all(reproducer.inspect_figure(figure_id) for figure_id in selected_ids)
            else:
                runs = render_figures(
                    figures,
                    selected_ids,
                    output_root,
                    jobs=args.jobs,
                    skip_existing=args.skip_existing,
                )
                errors = {run.figure_id: run.error for run in runs if run.error}
                for run in runs:
                    reproducer.merge_report(run.report)
                for figure_id, error in errors.items():
                    print(f"{figure_id}: {error}", file=sys.stderr)
                if args.check:
                    ok = check_outputs(figures, selected_ids, output_root, reproducer, errors)
                else:
                    ok = not errors and not reproducer.skipped_missing_inputs and not reproducer.failed
            reproducer.write_report(missing_report)
            return 0 if ok else 1
        finally:
            reproducer.close()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
