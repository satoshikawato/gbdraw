#!/usr/bin/env python3
"""Regenerate or verify every implemented documentation scenario."""

from __future__ import annotations

import argparse
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

from PIL import Image, ImageChops
from playwright.sync_api import BrowserType, sync_playwright


CAPTURE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = CAPTURE_ROOT.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(CAPTURE_ROOT) not in sys.path:
    sys.path.insert(0, str(CAPTURE_ROOT))

from config import (  # noqa: E402
    IMPLEMENTED_SCENARIO_IDS as GUI_SCENARIO_IDS,
    SUPPORTED_TIERS,
    first_circular_screenshot_paths,
    first_linear_screenshot_paths,
    gui_annotated_chloroplast_screenshot_paths,
    gui_annotation_tracks_screenshot_paths,
    gui_circular_layout_screenshot_paths,
    gui_circular_rings_screenshot_paths,
    gui_feature_highlight_screenshot_paths,
    gui_feature_presentation_screenshot_paths,
    gui_exports_screenshot_paths,
    gui_inputs_screenshot_paths,
    gui_interactive_editing_screenshot_paths,
    gui_interactive_handoff_screenshot_paths,
    gui_linear_layout_screenshot_paths,
    gui_losatn_screenshot_paths,
    gui_losatp_collinear_screenshot_paths,
    gui_losatp_tutorial_collinear_screenshot_paths,
    gui_losatp_groups_how_to_screenshot_paths,
    gui_losatp_groups_screenshot_paths,
    gui_quantitative_map_screenshot_paths,
    gui_quantitative_tracks_screenshot_paths,
    gui_precomputed_circular_rings_screenshot_paths,
    gui_session_reproduction_screenshot_paths,
    gui_styling_screenshot_paths,
    gui_tlosatx_screenshot_paths,
    gui_uploaded_comparison_screenshot_paths,
)
from docs.recipes.run_cli_scenarios import (  # noqa: E402
    IMPLEMENTED_SCENARIOS as CLI_SCENARIO_IDS,
    run_scenario as run_cli_scenario,
)
from docs.recipes.run_python_scenarios import (  # noqa: E402
    IMPLEMENTED_SCENARIOS as PYTHON_SCENARIO_IDS,
    run_scenario as run_python_scenario,
)
from flows.how_to.comparisons import (  # noqa: E402
    capture_gui_circular_rings,
    capture_gui_losatp_collinear,
    capture_gui_losatp_groups_how_to,
    capture_gui_tlosatx,
    capture_gui_uploaded_comparison,
)
from flows.how_to.inputs import capture_gui_inputs  # noqa: E402
from flows.how_to.exports import capture_gui_exports  # noqa: E402
from flows.how_to.interactive_sessions import (  # noqa: E402
    capture_gui_interactive_editing,
    capture_gui_session_reproduction,
)
from flows.how_to.layouts import (  # noqa: E402
    capture_gui_circular_layout,
    capture_gui_linear_layout,
)
from flows.how_to.presentation import (  # noqa: E402
    capture_gui_feature_presentation,
    capture_gui_styling,
)
from flows.how_to.tracks import (  # noqa: E402
    capture_gui_annotation_tracks,
    capture_gui_quantitative_tracks,
)
from flows.tutorials.gui_first_circular import capture_first_circular  # noqa: E402
from flows.tutorials.gui_first_linear import capture_first_linear  # noqa: E402
from flows.tutorials.gui_annotated_chloroplast import (  # noqa: E402
    capture_gui_annotated_chloroplast,
)
from flows.tutorials.gui_interactive_handoff import (  # noqa: E402
    capture_gui_interactive_handoff,
)
from flows.tutorials.gui_feature_highlight import (  # noqa: E402
    capture_gui_feature_highlight,
)
from flows.tutorials.gui_losatn import capture_gui_losatn  # noqa: E402
from flows.tutorials.gui_losatp_groups import (  # noqa: E402
    capture_gui_losatp_groups,
)
from flows.tutorials.gui_losatp_collinear import (  # noqa: E402
    capture_gui_losatp_tutorial_collinear,
)
from flows.tutorials.gui_precomputed_circular_rings import (  # noqa: E402
    capture_gui_precomputed_circular_rings,
)
from flows.tutorials.gui_quantitative_map import (  # noqa: E402
    capture_gui_quantitative_map,
)
from web_server import CaptureWebServer  # noqa: E402


CaptureFunction = Callable[
    [BrowserType, str, Mapping[str, Path], Path],
    Any,
]


@dataclass(frozen=True)
class ScenarioCapture:
    screenshot_paths: Callable[[], dict[str, Path]]
    capture: CaptureFunction
    tier: str


SCENARIOS = {
    "T-GUI-01": ScenarioCapture(
        screenshot_paths=first_circular_screenshot_paths,
        capture=capture_first_circular,
        tier="core",
    ),
    "T-GUI-02": ScenarioCapture(
        screenshot_paths=first_linear_screenshot_paths,
        capture=capture_first_linear,
        tier="core",
    ),
    "H-GUI-01": ScenarioCapture(
        screenshot_paths=gui_inputs_screenshot_paths,
        capture=capture_gui_inputs,
        tier="core",
    ),
    "T-GUI-03": ScenarioCapture(
        screenshot_paths=gui_losatn_screenshot_paths,
        capture=capture_gui_losatn,
        tier="extended",
    ),
    "H-GUI-02": ScenarioCapture(
        screenshot_paths=gui_circular_layout_screenshot_paths,
        capture=capture_gui_circular_layout,
        tier="extended",
    ),
    "H-GUI-03": ScenarioCapture(
        screenshot_paths=gui_linear_layout_screenshot_paths,
        capture=capture_gui_linear_layout,
        tier="extended",
    ),
    "H-GUI-04": ScenarioCapture(
        screenshot_paths=gui_uploaded_comparison_screenshot_paths,
        capture=capture_gui_uploaded_comparison,
        tier="extended",
    ),
    "H-GUI-05": ScenarioCapture(
        screenshot_paths=gui_tlosatx_screenshot_paths,
        capture=capture_gui_tlosatx,
        tier="extended",
    ),
    "H-GUI-06": ScenarioCapture(
        screenshot_paths=gui_circular_rings_screenshot_paths,
        capture=capture_gui_circular_rings,
        tier="extended",
    ),
    "T-GUI-04": ScenarioCapture(
        screenshot_paths=gui_losatp_groups_screenshot_paths,
        capture=capture_gui_losatp_groups,
        tier="extended",
    ),
    "T-GUI-05": ScenarioCapture(
        screenshot_paths=gui_annotated_chloroplast_screenshot_paths,
        capture=capture_gui_annotated_chloroplast,
        tier="extended",
    ),
    "T-GUI-06": ScenarioCapture(
        screenshot_paths=gui_precomputed_circular_rings_screenshot_paths,
        capture=capture_gui_precomputed_circular_rings,
        tier="extended",
    ),
    "T-GUI-08": ScenarioCapture(
        screenshot_paths=gui_losatp_tutorial_collinear_screenshot_paths,
        capture=capture_gui_losatp_tutorial_collinear,
        tier="extended",
    ),
    "T-GUI-09": ScenarioCapture(
        screenshot_paths=gui_interactive_handoff_screenshot_paths,
        capture=capture_gui_interactive_handoff,
        tier="extended",
    ),
    "T-GUI-10": ScenarioCapture(
        screenshot_paths=gui_feature_highlight_screenshot_paths,
        capture=capture_gui_feature_highlight,
        tier="extended",
    ),
    "T-GUI-12": ScenarioCapture(
        screenshot_paths=gui_quantitative_map_screenshot_paths,
        capture=capture_gui_quantitative_map,
        tier="extended",
    ),
    "H-GUI-07": ScenarioCapture(
        screenshot_paths=gui_losatp_groups_how_to_screenshot_paths,
        capture=capture_gui_losatp_groups_how_to,
        tier="extended",
    ),
    "H-GUI-08": ScenarioCapture(
        screenshot_paths=gui_losatp_collinear_screenshot_paths,
        capture=capture_gui_losatp_collinear,
        tier="extended",
    ),
    "H-GUI-09": ScenarioCapture(
        screenshot_paths=gui_quantitative_tracks_screenshot_paths,
        capture=capture_gui_quantitative_tracks,
        tier="extended",
    ),
    "H-GUI-10": ScenarioCapture(
        screenshot_paths=gui_annotation_tracks_screenshot_paths,
        capture=capture_gui_annotation_tracks,
        tier="extended",
    ),
    "H-GUI-11": ScenarioCapture(
        screenshot_paths=gui_styling_screenshot_paths,
        capture=capture_gui_styling,
        tier="extended",
    ),
    "H-GUI-12": ScenarioCapture(
        screenshot_paths=gui_feature_presentation_screenshot_paths,
        capture=capture_gui_feature_presentation,
        tier="extended",
    ),
    "H-GUI-13": ScenarioCapture(
        screenshot_paths=gui_interactive_editing_screenshot_paths,
        capture=capture_gui_interactive_editing,
        tier="extended",
    ),
    "H-GUI-14": ScenarioCapture(
        screenshot_paths=gui_session_reproduction_screenshot_paths,
        capture=capture_gui_session_reproduction,
        tier="extended",
    ),
    "H-GUI-15": ScenarioCapture(
        screenshot_paths=gui_exports_screenshot_paths,
        capture=capture_gui_exports,
        tier="core",
    ),
}

TIER_RANK = {tier: index for index, tier in enumerate(SUPPORTED_TIERS)}
ALL_SCENARIO_IDS = (*GUI_SCENARIO_IDS, *CLI_SCENARIO_IDS, *PYTHON_SCENARIO_IDS)
MAX_RASTER_NOISE_PIXELS = 100
MAX_RASTER_CHANNEL_DELTA = 1

if len(ALL_SCENARIO_IDS) != len(set(ALL_SCENARIO_IDS)):
    raise RuntimeError("Documentation scenario IDs must be unique across surfaces.")


def _images_match(expected_path: Path, actual_path: Path) -> bool:
    with Image.open(expected_path) as expected, Image.open(actual_path) as actual:
        if expected.size != actual.size or expected.mode != actual.mode:
            return False
        difference = ImageChops.difference(expected, actual)
        if difference.getbbox() is None:
            return True

        extrema = difference.getextrema()
        if any(high > MAX_RASTER_CHANNEL_DELTA for _, high in extrema):
            return False
        changed_pixels = sum(
            pixel != (0, 0, 0) for pixel in difference.getdata()
        )
        return changed_pixels <= MAX_RASTER_NOISE_PIXELS


def _capture_scenario(
    scenario_id: str,
    output_paths: Mapping[str, Path],
    browser_type: BrowserType,
    base_url: str,
    download_root: Path,
) -> None:
    result = SCENARIOS[scenario_id].capture(
        browser_type,
        base_url,
        output_paths,
        download_root / scenario_id.lower(),
    )
    total_screenshot_bytes = sum(result.screenshot_bytes.values())
    print(
        f"{scenario_id}: captured {len(result.screenshot_bytes)} images "
        f"({total_screenshot_bytes} bytes, "
        f"{result.final_svg_semantics['featureElementCount']} features)"
    )
    print(
        f"{scenario_id}: validated "
        f"{result.download['filename']} ({result.download['bytes']} bytes)"
    )
    if getattr(result, "tsv_download", None):
        print(
            f"{scenario_id}: validated "
            f"{result.tsv_download['filename']} "
            f"({result.tsv_download['bytes']} bytes, "
            f"{result.tsv_download['rows']} rows)"
        )
    if getattr(result, "fasta_download", None):
        print(
            f"{scenario_id}: validated "
            f"{result.fasta_download['filename']} "
            f"({result.fasta_download['bytes']} bytes, "
            f"{result.fasta_download['records']} records)"
        )
        verified_members = result.fasta_download.get("verifiedMembers", [])
        if verified_members:
            summary = ", ".join(
                f"{member['proteinId']}={member['length']}aa:"
                f"{member['sha256'][:12]}"
                for member in verified_members
            )
            print(
                f"{scenario_id}: FASTA matches raw GenBank CDS translations "
                f"({summary})"
            )
    if getattr(result, "session_download", None):
        print(
            f"{scenario_id}: validated "
            f"{result.session_download['filename']} "
            f"({result.session_download['bytes']} bytes, "
            f"session v{result.session_download['version']})"
        )
    if getattr(result, "downloads", None):
        for artifact_name, artifact in result.downloads.items():
            if artifact.get("filename") == result.download.get("filename"):
                continue
            print(
                f"{scenario_id}: validated {artifact_name} "
                f"{artifact['filename']} ({artifact['bytes']} bytes)"
            )


def _capture_many(
    scenario_ids: tuple[str, ...],
    output_paths_by_scenario: Mapping[str, Mapping[str, Path]],
) -> None:
    with (
        tempfile.TemporaryDirectory(prefix="gbdraw-doc-download-") as download_dir,
        CaptureWebServer() as server,
        sync_playwright() as playwright,
    ):
        download_root = Path(download_dir)
        for scenario_id in scenario_ids:
            _capture_scenario(
                scenario_id,
                output_paths_by_scenario[scenario_id],
                playwright.chromium,
                server.base_url,
                download_root,
            )


def _check(
    scenario_ids: tuple[str, ...],
    committed_paths_by_scenario: Mapping[str, Mapping[str, Path]],
) -> None:
    missing = [
        str(path)
        for scenario_id in scenario_ids
        for path in committed_paths_by_scenario[scenario_id].values()
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            f"Missing committed screenshots for --check: {', '.join(missing)}"
        )

    with tempfile.TemporaryDirectory(prefix="gbdraw-doc-capture-") as temp_dir:
        candidate_paths_by_scenario = {
            scenario_id: {
                name: Path(temp_dir) / scenario_id.lower() / name
                for name in committed_paths_by_scenario[scenario_id]
            }
            for scenario_id in scenario_ids
        }
        _capture_many(scenario_ids, candidate_paths_by_scenario)

        for scenario_id in scenario_ids:
            committed_paths = committed_paths_by_scenario[scenario_id]
            candidate_paths = candidate_paths_by_scenario[scenario_id]
            stale = [
                name
                for name, committed_path in committed_paths.items()
                if not _images_match(committed_path, candidate_paths[name])
            ]
            if stale:
                raise AssertionError(
                    f"The committed {scenario_id} screenshots are stale: "
                    f"{', '.join(stale)}. Regenerate with "
                    "python docs/capture/run_all.py "
                    f"--scenario {scenario_id} --tier {SCENARIOS[scenario_id].tier}"
                )
            print(f"{scenario_id}: all committed screenshots match a fresh capture")


def _run_recipes(
    scenario_ids: tuple[str, ...],
    run_scenario: Callable[..., tuple[Path, ...]],
    *,
    check: bool,
) -> None:
    action = "verified" if check else "wrote"
    for scenario_id in scenario_ids:
        for destination in run_scenario(scenario_id, check=check):
            print(f"{scenario_id}: {action} {destination}")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate deterministic gbdraw documentation artifacts."
    )
    parser.add_argument(
        "--scenario",
        choices=("all", *ALL_SCENARIO_IDS),
        default="all",
        help="Run one GUI, CLI, or Python scenario, or every implemented scenario.",
    )
    parser.add_argument(
        "--tier",
        choices=SUPPORTED_TIERS,
        default="core",
        help=(
            "GUI capture tier. Higher tiers include lower-tier GUI scenarios; "
            "all CLI and Python recipes always run when --scenario is all."
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Regenerate each selected artifact and compare it with the committed output.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    included_gui_ids = tuple(
        scenario_id
        for scenario_id in GUI_SCENARIO_IDS
        if TIER_RANK[SCENARIOS[scenario_id].tier] <= TIER_RANK[args.tier]
    )

    if args.scenario == "all" or args.scenario in SCENARIOS:
        if args.scenario != "all" and args.scenario not in included_gui_ids:
            raise ValueError(
                f"{args.scenario} requires --tier {SCENARIOS[args.scenario].tier}"
            )
        gui_ids = (
            included_gui_ids if args.scenario == "all" else (args.scenario,)
        )
        output_paths_by_scenario = {
            scenario_id: SCENARIOS[scenario_id].screenshot_paths()
            for scenario_id in gui_ids
        }
        if args.check:
            _check(gui_ids, output_paths_by_scenario)
        else:
            _capture_many(gui_ids, output_paths_by_scenario)

    if args.scenario == "all" or args.scenario in CLI_SCENARIO_IDS:
        cli_ids = CLI_SCENARIO_IDS if args.scenario == "all" else (args.scenario,)
        _run_recipes(cli_ids, run_cli_scenario, check=args.check)

    if args.scenario == "all" or args.scenario in PYTHON_SCENARIO_IDS:
        python_ids = (
            PYTHON_SCENARIO_IDS if args.scenario == "all" else (args.scenario,)
        )
        _run_recipes(python_ids, run_python_scenario, check=args.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
