#!/usr/bin/env python
"""Capture compact circular preset radial geometry fixtures.

This script records resolver geometry only. It intentionally avoids full SVG
output and font-dependent external label measurements.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from Bio import SeqIO

from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import DepthConfigurator
from gbdraw.diagrams.circular.presets import CircularPresetContext, circular_radial_plan_for_preset
from gbdraw.diagrams.circular.radial_layout import RadialBand, resolve_circular_radial_layout
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.features.shapes import resolve_directional_feature_types
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering


BASE_COMMIT = "5dd1b69fe066bc0b39cff7e072bec003c6e654ae"
DEFAULT_INPUTS = (
    "tests/test_inputs/MG1655.gbk",
    "tests/test_inputs/AP027280.gb",
)
SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]
VISIBILITY_MATRIX = {
    "features_ticks": {"show_depth": False, "show_gc": False, "show_skew": False},
    "gc_skew": {"show_depth": False, "show_gc": True, "show_skew": True},
    "depth_gc_skew": {"show_depth": True, "show_gc": True, "show_skew": True},
}


def _band_data(band: RadialBand | None) -> dict[str, float] | None:
    if band is None:
        return None
    return {
        "inner_px": float(band.inner_px),
        "outer_px": float(band.outer_px),
        "width_px": float(band.width_px),
        "center_px": float(band.center_px),
    }


def _slot_data(slot: Any) -> dict[str, Any]:
    return {
        "slot_index": int(slot.slot_index),
        "id": str(slot.id),
        "renderer": str(slot.renderer),
        "side": str(slot.side),
        "z": int(slot.z),
        "anchor_radius_px": None if slot.anchor_radius_px is None else float(slot.anchor_radius_px),
        "anchor_offset_px": None if slot.anchor_offset_px is None else float(slot.anchor_offset_px),
        "requested_width_px": None if slot.requested_width_px is None else float(slot.requested_width_px),
        "resolved_width_px": None if slot.resolved_width_px is None else float(slot.resolved_width_px),
        "packing_band_px": _band_data(slot.packing_band_px),
        "draw_band_px": _band_data(slot.draw_band_px),
        "reserved_band_px": _band_data(slot.reserved_band_px),
        "explicit_anchor": bool(slot.explicit_anchor),
        "explicit_width": bool(slot.explicit_width),
        "compressed": bool(slot.compressed),
        "params": dict(slot.params or {}),
    }


def _feature_dict(record: Any, cfg: GbdrawConfig) -> dict[str, Any]:
    color_table, default_colors = preprocess_color_tables(None, load_default_colors("", "default"))
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    feature_dict, _used = create_feature_dict(
        record,
        color_table,
        SELECTED_FEATURES,
        default_colors,
        cfg.canvas.strandedness,
        cfg.canvas.resolve_overlaps,
        label_filtering,
        directional_feature_types=resolve_directional_feature_types(None),
        feature_visibility_rules=None,
        compute_label_text=False,
    )
    return feature_dict


def capture_case(
    *,
    input_path: Path,
    preset: str,
    strandedness: bool,
    visibility_name: str,
    visibility: dict[str, bool],
) -> dict[str, Any]:
    record = SeqIO.read(str(input_path), "genbank")
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=visibility["show_gc"],
        show_skew=visibility["show_skew"],
        show_depth=visibility["show_depth"],
        track_type=preset,
        strandedness=strandedness,
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = CircularCanvasConfigurator(
        input_path.stem,
        config_dict,
        "none",
        record,
        cfg=cfg,
    )
    context = CircularPresetContext(
        cfg=cfg,
        canvas_config=canvas_config,
        total_length=len(record.seq),
        strandedness=strandedness,
        show_features=True,
        show_ticks=True,
        show_depth=visibility["show_depth"],
        show_gc=visibility["show_gc"],
        show_skew=visibility["show_skew"],
        dinucleotide="GC",
    )
    radial_plan = circular_radial_plan_for_preset(preset, context)
    depth_config = (
        DepthConfigurator(
            1,
            1,
            config_dict,
            cfg=cfg,
        )
        if visibility["show_depth"]
        else None
    )
    layout = resolve_circular_radial_layout(
        total_length=len(record.seq),
        canvas_config=canvas_config,
        cfg=cfg,
        slots=radial_plan.slots,
        feature_dict=_feature_dict(record, cfg),
        show_features=True,
        show_ticks=True,
        preferred_anchor_slot_ids=radial_plan.preferred_anchor_slot_ids,
        depth_config=depth_config,
    )
    return {
        "input": str(input_path),
        "preset": preset,
        "strandedness": strandedness,
        "visibility": visibility_name,
        "axis_radius_px": float(layout.axis.radius_px),
        "outer_content_radius_px": float(layout.outer_content_radius_px),
        "slots": [_slot_data(slot) for slot in layout.slots],
    }


def capture(inputs: list[Path]) -> dict[str, Any]:
    cases: list[dict[str, Any]] = []
    for input_path in inputs:
        for preset in ("tuckin", "middle", "spreadout"):
            for strandedness in (True, False):
                for visibility_name, visibility in VISIBILITY_MATRIX.items():
                    cases.append(
                        capture_case(
                            input_path=input_path,
                            preset=preset,
                            strandedness=strandedness,
                            visibility_name=visibility_name,
                            visibility=visibility,
                        )
                    )
    return {
        "base_commit": BASE_COMMIT,
        "capture_source": "current_worktree",
        "schema": 1,
        "description": (
            "Compact circular preset radial geometry oracle. Run this script from "
            "the target checkout when regenerating expected geometry."
        ),
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/circular_preset_oracle") / f"{BASE_COMMIT}.json",
        help="Output JSON fixture path.",
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        default=[Path(path) for path in DEFAULT_INPUTS],
        help="GenBank inputs to capture.",
    )
    args = parser.parse_args()

    data = capture(list(args.inputs))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
