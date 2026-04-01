from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import Any

from Bio import SeqIO
from PIL import Image

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.reproduce_examples_manifest import (  # noqa: E402
    BlastPreparation,
    CliRecipe,
    CompositePanel,
    CompositeRecipe,
    FigureSpec,
    FastaPreparation,
    FileArgument,
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
        default=PROJECT_ROOT / "_reproduced",
        help="Root directory for reproduced outputs (default: %(default)s).",
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
        self.examples_dir = project_root / "examples"
        self.test_inputs_dir = project_root / "tests" / "test_inputs"
        self.output_root = output_root
        self.support_dir = output_root / "_support_assets"
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
        self._prepared_cache: dict[str, Path] = {}
        self._active_render: set[str] = set()
        self._active_inspect: set[str] = set()
        self._tmpdir = tempfile.TemporaryDirectory(prefix="reproduce_examples_")
        self.temp_root = Path(self._tmpdir.name)

    def close(self) -> None:
        self._tmpdir.cleanup()

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
                    "could_derive_if_base_inputs_present": data["could_derive_if_base_inputs_present"],
                    "base_inputs": sorted(data["base_inputs"]),
                }
            )
        return {
            "generated": self.generated,
            "skipped_missing_inputs": self.skipped_missing_inputs,
            "failed": self.failed,
            "aliases_used": alias_entries,
            "missing_inputs": missing_entries,
        }

    def write_report(self, report_path: Path) -> None:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(self.report_payload(), indent=2, sort_keys=False) + "\n")

    def add_missing(
        self,
        filename: str,
        figure_id: str,
        *,
        could_derive: bool,
        base_inputs: tuple[str, ...] = (),
    ) -> None:
        entry = self.missing_inputs.setdefault(
            filename,
            {
                "figures": set(),
                "could_derive_if_base_inputs_present": False,
                "base_inputs": set(),
            },
        )
        entry["figures"].add(figure_id)
        if could_derive:
            entry["could_derive_if_base_inputs_present"] = True
        entry["base_inputs"].update(base_inputs)

    def add_alias_use(self, requested: str, resolved: Path, figure_id: str) -> None:
        resolved_str = str(resolved.relative_to(self.project_root))
        self.aliases_used[(requested, resolved_str)].add(figure_id)

    def _existing_input_path(self, filename: str) -> Path | None:
        for directory in (self.examples_dir, self.test_inputs_dir, self.support_dir):
            candidate = directory / filename
            if candidate.exists():
                return candidate
        return None

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
        if not target_path.exists():
            target_path.write_text(asset.content)
        return target_path

    def _base_inputs_for(
        self,
        filename: str,
        preparation_map: dict[str, FastaPreparation | BlastPreparation],
        seen: set[str] | None = None,
    ) -> tuple[str, ...]:
        if seen is None:
            seen = set()
        if filename in seen:
            return ()
        seen.add(filename)
        preparation = preparation_map.get(filename)
        if preparation is None:
            return (filename,)
        if isinstance(preparation, FastaPreparation):
            return self._base_inputs_for(preparation.source_filename, preparation_map, seen)
        query_inputs = self._base_inputs_for(preparation.query_filename, preparation_map, seen)
        subject_inputs = self._base_inputs_for(preparation.subject_filename, preparation_map, seen)
        return tuple(dict.fromkeys(query_inputs + subject_inputs))

    def resolve_input(
        self,
        filename: str,
        figure_id: str,
        preparation_map: dict[str, FastaPreparation | BlastPreparation],
        *,
        dry_run: bool,
        stack: tuple[str, ...] = (),
    ) -> Path | None:
        cached_path = self._prepared_cache.get(filename)
        if cached_path and cached_path.exists():
            return cached_path

        existing = self._existing_input_path(filename)
        if existing:
            return existing

        alias = self._alias_input_path(filename, figure_id)
        if alias:
            return alias

        asset_path = self._materialize_support_asset(filename, dry_run)
        if asset_path:
            return asset_path

        if filename in preparation_map:
            return self.ensure_preparation(
                preparation_map[filename],
                figure_id,
                preparation_map,
                dry_run=dry_run,
                stack=stack,
            )

        self.add_missing(filename, figure_id, could_derive=False)
        return None

    def ensure_preparation(
        self,
        preparation: FastaPreparation | BlastPreparation,
        figure_id: str,
        preparation_map: dict[str, FastaPreparation | BlastPreparation],
        *,
        dry_run: bool,
        stack: tuple[str, ...] = (),
    ) -> Path | None:
        if preparation.output_filename in stack:
            raise RuntimeError(f"Cyclic preparation detected for {preparation.output_filename}")

        target_path = self.support_dir / preparation.output_filename
        if target_path.exists():
            self._prepared_cache[preparation.output_filename] = target_path
            return target_path

        next_stack = stack + (preparation.output_filename,)
        if isinstance(preparation, FastaPreparation):
            source_path = self.resolve_input(
                preparation.source_filename,
                figure_id,
                preparation_map,
                dry_run=dry_run,
                stack=next_stack,
            )
            if source_path is None:
                self.add_missing(
                    preparation.output_filename,
                    figure_id,
                    could_derive=True,
                    base_inputs=(preparation.source_filename,),
                )
                return None
            if dry_run:
                return target_path
            target_path.parent.mkdir(parents=True, exist_ok=True)
            records = list(SeqIO.parse(str(source_path), "genbank"))
            if not records:
                raise RuntimeError(f"No GenBank records found in {source_path}")
            with target_path.open("w", encoding="utf-8") as handle:
                SeqIO.write(records, handle, "fasta")
            self._prepared_cache[preparation.output_filename] = target_path
            return target_path

        query_path = self.resolve_input(
            preparation.query_filename,
            figure_id,
            preparation_map,
            dry_run=dry_run,
            stack=next_stack,
        )
        subject_path = self.resolve_input(
            preparation.subject_filename,
            figure_id,
            preparation_map,
            dry_run=dry_run,
            stack=next_stack,
        )
        if query_path is None or subject_path is None:
            self.add_missing(
                preparation.output_filename,
                figure_id,
                could_derive=True,
                base_inputs=self._base_inputs_for(preparation.output_filename, preparation_map),
            )
            return None
        if dry_run:
            if shutil.which(preparation.tool) is None:
                return None
            return target_path
        if shutil.which(preparation.tool) is None:
            raise RuntimeError(f"Required BLAST tool not found on PATH: {preparation.tool}")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            preparation.tool,
            "-query",
            str(query_path),
            "-subject",
            str(subject_path),
            "-outfmt",
            preparation.outfmt,
            "-out",
            str(target_path),
        ]
        result = subprocess.run(
            cmd,
            cwd=str(self.project_root),
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            raise RuntimeError(f"{preparation.tool} failed for {preparation.output_filename}: {result.stderr.strip()}")
        self._prepared_cache[preparation.output_filename] = target_path
        return target_path

    def _resolve_cli_recipe(
        self,
        recipe: CliRecipe,
        figure_id: str,
        preparation_map: dict[str, FastaPreparation | BlastPreparation],
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
                path = self.resolve_input(filename, figure_id, preparation_map, dry_run=dry_run)
                if path is None:
                    return None
                resolved[key].append(path)
        for file_arg in recipe.file_args:
            resolved_paths: list[Path] = []
            for filename in file_arg.filenames:
                path = self.resolve_input(filename, figure_id, preparation_map, dry_run=dry_run)
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
        cmd = [sys.executable, "-m", "gbdraw.cli", recipe.subcommand]
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
        cmd.extend(["-o", str(output_prefix), "-f", fmt])
        cmd.extend(resolved["extra_args"])
        return cmd

    def _run_cli_recipe(
        self,
        recipe: CliRecipe,
        figure_id: str,
        preparation_map: dict[str, FastaPreparation | BlastPreparation],
        target_path: Path,
        *,
        dry_run: bool,
    ) -> bool:
        resolved = self._resolve_cli_recipe(recipe, figure_id, preparation_map, dry_run=dry_run)
        if resolved is None:
            return False
        if dry_run:
            return True

        target_path.parent.mkdir(parents=True, exist_ok=True)
        output_prefix = target_path.with_suffix("")
        cmd = self._build_command(recipe, resolved, output_prefix, target_path.suffix.lstrip("."))
        result = subprocess.run(
            cmd,
            cwd=str(self.project_root),
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or result.stdout.strip() or f"gbdraw {recipe.subcommand} failed")
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

    def _load_image(self, source_path: Path, target_size: tuple[int, int]) -> Image.Image:
        if source_path.suffix.lower() == ".svg":
            try:
                import cairosvg
            except ModuleNotFoundError as exc:  # pragma: no cover - environment dependent
                raise RuntimeError("CairoSVG is required to build composite PNG outputs") from exc
            src_width, src_height = self._svg_size(source_path)
            scale = min(target_size[0] / src_width, target_size[1] / src_height)
            render_width = max(1, int(round(src_width * scale)))
            render_height = max(1, int(round(src_height * scale)))
            png_bytes = cairosvg.svg2png(
                url=str(source_path),
                output_width=render_width,
                output_height=render_height,
            )
            image = Image.open(BytesIO(png_bytes)).convert("RGBA")
        else:
            image = Image.open(source_path).convert("RGBA")
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
        *,
        dry_run: bool,
    ) -> Image.Image | None:
        if panel.figure_id:
            if dry_run:
                return None if not self.inspect_figure(panel.figure_id) else Image.new("RGBA", target_size, (255, 255, 255, 0))
            if not self.render_figure(panel.figure_id):
                return None
            return self._load_image(self.output_path_for(panel.figure_id), target_size)
        if panel.recipe is None:
            raise RuntimeError(f"Composite panel {panel_index} in {figure_id} has no source")
        prep_map = {prep.output_filename: prep for prep in panel.preparations}
        if dry_run:
            return None if not self._run_cli_recipe(panel.recipe, figure_id, prep_map, self.temp_root / "dry_run.svg", dry_run=True) else Image.new("RGBA", target_size, (255, 255, 255, 0))
        panel_dir = self.temp_root / figure_id
        panel_dir.mkdir(parents=True, exist_ok=True)
        panel_path = panel_dir / f"panel_{panel_index}.svg"
        if not self._run_cli_recipe(panel.recipe, figure_id, prep_map, panel_path, dry_run=False):
            return None
        return self._load_image(panel_path, target_size)

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
                image = self._render_panel_source(panel, figure_id, index, recipe.tile_size, dry_run=True)
                if image is None:
                    return False
            return True

        panel_images: list[Image.Image] = []
        for index, panel in enumerate(recipe.panels):
            image = self._render_panel_source(panel, figure_id, index, recipe.tile_size, dry_run=False)
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
        target_path.parent.mkdir(parents=True, exist_ok=True)
        background.convert("RGB").save(target_path)
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
        for index, panel in enumerate(recipe.panels):
            if panel.figure_id is None or panel.box is None:
                raise RuntimeError(f"social_preview panel {index} for {figure_id} must reference a figure and a box")
            if not self.render_figure(panel.figure_id):
                return False
            box_width = panel.box[2]
            box_height = panel.box[3]
            image = self._load_image(self.output_path_for(panel.figure_id), (box_width, box_height))
            canvas.paste(image, (panel.box[0], panel.box[1]), image)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        canvas.convert("RGB").save(target_path)
        return True

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
                prep_map = {prep.output_filename: prep for prep in spec.preparations}
                ready = self._run_cli_recipe(spec.recipe, figure_id, prep_map, target_path, dry_run=True)
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
                prep_map = {prep.output_filename: prep for prep in spec.preparations}
                success = self._run_cli_recipe(spec.recipe, figure_id, prep_map, target_path, dry_run=False)
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


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    figures = build_figure_specs()
    selected_ids = select_figures(figures, args.group, args.figure)
    output_root = args.output_root if args.output_root.is_absolute() else PROJECT_ROOT / args.output_root
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
            exit_code = 0 if all(reproducer.inspect_figure(figure_id) for figure_id in selected_ids) else 1
        else:
            exit_code = 0
            for figure_id in selected_ids:
                try:
                    if not reproducer.render_figure(figure_id):
                        exit_code = 1
                except Exception as exc:
                    print(f"{figure_id}: {exc}", file=sys.stderr)
                    exit_code = 1
        if reproducer.skipped_missing_inputs or reproducer.failed:
            exit_code = 1
        reproducer.write_report(missing_report)
        return exit_code
    finally:
        reproducer.close()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
