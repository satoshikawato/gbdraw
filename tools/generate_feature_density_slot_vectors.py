#!/usr/bin/env python3
"""Generate the shared feature density row vectors (Python's stacks, exact).

``tests/fixtures/feature_density_slot_vectors.json`` holds, per case, the
feature density track ids, the stack they are added to (one of the default
stacks, or the listed rows), and the rows Python gives it:
``gbdraw.tracks.parsing.with_feature_density_slots`` over
``default_circular_track_slots`` / ``default_linear_track_slots`` with the slot
the drawing assemblers make for each track.

- ``tests/test_feature_density_slot_vectors.py`` regenerates the file in memory
  (this generator's ``--check``).
- ``tests/web/feature-density-slot-vectors.test.mjs`` runs the Web's
  ``insertFeatureDensitySlots``, ``createDefaultCircularTrackSlots`` and
  ``createDefaultLinearTrackSlots`` on the same inputs and compares.

No case holds an Axis index: Python adds these rows only to a stack the reader did
not write, and an Axis index requires a written stack (``linear_track_axis_index
requires linear_track_slots`` in ``gbdraw/api/diagram.py``; Circular has none), so
the two never meet. The Web's own Axis shift is tested without Python vectors.

A row is compared by ``id``, ``renderer`` and ``track`` (the ``track`` param of
a feature density row; ``side`` too on a feature density row, null when unset).

Usage:
  python tools/generate_feature_density_slot_vectors.py            # rewrite the vectors
  python tools/generate_feature_density_slot_vectors.py --check    # fail when stale
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "tests" / "fixtures" / "feature_density_slot_vectors.json"
sys.path.insert(0, str(REPO_ROOT))

from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot  # noqa: E402
from gbdraw.tracks.circular import default_circular_track_slots  # noqa: E402
from gbdraw.tracks.linear import default_linear_track_slots  # noqa: E402
from gbdraw.tracks.parsing import feature_density_slot_id, with_feature_density_slots  # noqa: E402

TRACK_SETS: tuple[list[str], ...] = ([], ["genes"], ["genes", "repeats"])
STACK_ROWS: dict[str, list[tuple[str, str]]] = {
    "features_first": [("features", "features"), ("gc_content", "dinucleotide_content"), ("gc_skew", "dinucleotide_skew")],
    "no_features": [("gc_content", "dinucleotide_content"), ("gc_skew", "dinucleotide_skew")],
    "features_last": [("gc_content", "dinucleotide_content"), ("features", "features")],
    "already_drawn": [("features", "features"), ("density_old", "feature_density"), ("gc_content", "dinucleotide_content")],
}
CIRCULAR_DEFAULTS: tuple[dict[str, Any], ...] = (
    {},
    {"showDepth": True, "depthTrackCount": 2},
    {"showTicks": False, "showGc": False, "showSkew": False},
)
LINEAR_DEFAULTS: tuple[dict[str, Any], ...] = (
    {},
    {"showGc": True, "showSkew": True},
    {"showDepth": True, "depthTrackCount": 2, "showGc": True, "trackLayout": "above"},
)


def _slot_class(mode: str) -> type:
    return CircularTrackSlot if mode == "circular" else LinearTrackSlot


def _density_side(mode: str) -> dict[str, str]:
    """Linear assembly puts the row below the Axis; Circular leaves the side unset."""

    return {"side": "below"} if mode == "linear" else {}


def _make_density_slot(mode: str, track_id: str) -> Any:
    return _slot_class(mode)(
        id=feature_density_slot_id(track_id), renderer="feature_density", params={"track": track_id}, **_density_side(mode)
    )


def _given_slot(mode: str, row_id: str, renderer: str) -> Any:
    if renderer == "feature_density":
        return _slot_class(mode)(id=row_id, renderer=renderer, params={"track": "old"}, **_density_side(mode))
    return _slot_class(mode)(id=row_id, renderer=renderer)


def _row(slot: Any) -> dict[str, Any]:
    row: dict[str, Any] = {"id": str(slot.id), "renderer": str(slot.renderer)}
    if slot.renderer == "feature_density":
        row["track"] = str(slot.params["track"])
        row["side"] = slot.side
    return row


def _default_stack(mode: str, options: dict[str, Any]) -> list[Any]:
    depth = {"show_depth": options.get("showDepth", False), "depth_track_count": options.get("depthTrackCount", 1)}
    if mode == "circular":
        return default_circular_track_slots(
            show_ticks=options.get("showTicks", True), show_gc=options.get("showGc", True),
            show_skew=options.get("showSkew", True), **depth,
        )
    return default_linear_track_slots(
        show_gc=options.get("showGc", False), show_skew=options.get("showSkew", False),
        track_layout=options.get("trackLayout", "middle"), **depth,
    )


def _case_name(mode: str, stack: str, tracks: list[str]) -> str:
    return f"{mode}-{stack}-{len(tracks)}tracks"


def cases() -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for mode, defaults in (("circular", CIRCULAR_DEFAULTS), ("linear", LINEAR_DEFAULTS)):
        for tracks in TRACK_SETS:
            for number, options in enumerate(defaults):
                slots = with_feature_density_slots(
                    _default_stack(mode, options), tracks, lambda track_id, mode=mode: _make_density_slot(mode, track_id)
                )
                out.append({
                    "id": _case_name(mode, f"default{number}", tracks), "mode": mode, "kind": "default",
                    "options": options, "trackIds": tracks, "expected": {"rows": [_row(slot) for slot in slots]},
                })
            for stack, rows in STACK_ROWS.items():
                given = [_given_slot(mode, row_id, renderer) for row_id, renderer in rows]
                slots = with_feature_density_slots(given, tracks, lambda track_id, mode=mode: _make_density_slot(mode, track_id))
                out.append({
                    "id": _case_name(mode, stack, tracks), "mode": mode, "kind": "stack",
                    "rows": [_row(slot) for slot in given], "trackIds": tracks,
                    "expected": {"rows": [_row(slot) for slot in slots]},
                })
    return out


def render() -> str:
    return json.dumps({"schemaVersion": 1, "cases": cases()}, indent=2, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Fail if the vectors are stale.")
    args = parser.parse_args()
    expected = render()
    if args.check:
        if not TARGET.is_file() or TARGET.read_text(encoding="utf-8") != expected:
            print("Feature density slot vectors are stale. Run: python tools/generate_feature_density_slot_vectors.py")
            return 1
        return 0
    TARGET.write_text(expected, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
