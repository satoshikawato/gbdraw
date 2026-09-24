"""Validate the fixed video timeline and captured source bundle."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


EXPECTED_SCENES = (
    "intro", "circular", "plastome", "genome-comparison", "cluster-comparison",
    "labels-product", "labels-gene", "colors-before", "colors-functional",
    "dloop-before", "dloop-bracket", "export-svg", "outro",
)
EXPECTED_ASSETS = frozenset({
    "human.circular", "tobacco.plastome", "lambda-de3.comparison",
    "bgc.comparison", "human.labels-product", "human.labels-gene",
    "human.functional-colors", "human.dloop-bracket", "human.svg-export",
})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contained_file(root: Path, relative: str) -> Path:
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts or not candidate.parts:
        raise ValueError(f"Unsafe relative path: {relative!r}")
    resolved = (root / candidate).resolve()
    if not resolved.is_relative_to(root.resolve()) or not resolved.is_file():
        raise ValueError(f"Missing or escaped file: {relative!r}")
    return resolved


def load_storyboard(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if (data.get("schema_version"), data.get("video_id"), data.get("width"),
            data.get("height"), data.get("fps"), data.get("total_frames")) != (
            1, "meet-gbdraw", 1920, 1080, 30, 1140):
        raise ValueError("Unsupported Meet gbdraw video format")
    scenes = data.get("scenes")
    if not isinstance(scenes, list) or tuple(row.get("id") for row in scenes) != EXPECTED_SCENES:
        raise ValueError("Storyboard must contain the 13 scenes in the approved order")
    total = 0
    for scene in scenes:
        assets = scene.get("assets")
        if not isinstance(scene.get("frames"), int) or scene["frames"] <= 0:
            raise ValueError(f"Invalid frame count for {scene['id']}")
        if not isinstance(assets, list) or not assets or any(a not in EXPECTED_ASSETS for a in assets):
            raise ValueError(f"Unknown or missing asset in {scene['id']}")
        if scene.get("layout") not in {"single", "grid2x2"} or len(assets) != (4 if scene["layout"] == "grid2x2" else 1):
            raise ValueError(f"Invalid scene layout for {scene['id']}")
        if not isinstance(scene.get("caption"), str) or not scene["caption"].strip():
            raise ValueError(f"Missing caption for {scene['id']}")
        scene["start_frame"] = total
        total += scene["frames"]
    if total != 1140:
        raise ValueError(f"Expected 1,140 frames, found {total}")
    return data


def load_assets(path: Path, *, complete: bool = True) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or not isinstance(data.get("assets"), dict):
        raise ValueError("Invalid assets manifest")
    if complete and set(data["assets"]) != EXPECTED_ASSETS:
        raise ValueError("The nine required source assets are incomplete")
    for asset_id, asset in data["assets"].items():
        if asset_id not in EXPECTED_ASSETS or asset.get("kind") not in {"image", "video"}:
            raise ValueError(f"Invalid asset {asset_id}")
        source = contained_file(path.parent, asset["path"])
        if sha256(source) != asset.get("sha256"):
            raise ValueError(f"Asset checksum mismatch: {asset_id}")
        for provenance in ("svg", "evidence"):
            if asset.get(provenance):
                original = contained_file(path.parent, asset[provenance]["path"])
                if sha256(original) != asset[provenance].get("sha256"):
                    raise ValueError(f"{provenance} checksum mismatch: {asset_id}")
    return data
