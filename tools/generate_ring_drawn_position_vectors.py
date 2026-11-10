#!/usr/bin/env python3
"""Generate the shared ring drawn-position vectors (Python's drawn order, exact).

``tests/fixtures/ring_drawn_position_vectors.json`` holds, per case, a number of
comparison rings, a track stack (or none), and the series places of the rings
Python draws, outermost first: ``_drawn_ring_colors`` in
``gbdraw/analysis/conservation.py``, the order the scheme colours follow
(PD-OI-110). A case without a stack draws over the slots Python makes for the
rings (``_conservation_slots_for_tracks`` inserted into the default stack by
``_insert_conservation_slots``).

- ``tests/test_ring_drawn_position_vectors.py`` regenerates the file in memory
  (this generator's ``--check``).
- ``tests/web/ring-drawn-position-vectors.test.mjs`` runs the Web's
  ``ringDrawnPositions`` on the same inputs: each drawn ring's position is its
  place in Python's list.

Usage:
  python tools/generate_ring_drawn_position_vectors.py            # rewrite the vectors
  python tools/generate_ring_drawn_position_vectors.py --check    # fail when stale
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "tests" / "fixtures" / "ring_drawn_position_vectors.json"
sys.path.insert(0, str(REPO_ROOT))

from gbdraw.analysis.conservation import ConservationTrack, _drawn_ring_colors  # noqa: E402
from gbdraw.api.diagram import _conservation_slots_for_tracks, _insert_conservation_slots  # noqa: E402
from gbdraw.config.models.root import GbdrawConfig  # noqa: E402
from gbdraw.config.toml import load_config_toml  # noqa: E402
from gbdraw.tracks import CircularTrackSlot  # noqa: E402
from gbdraw.tracks.circular import default_circular_track_slots  # noqa: E402


def _ring(track_index: Any, *, enabled: bool = True, slot_id: str | None = None) -> dict[str, Any]:
    return {"id": slot_id or f"ring_{track_index}", "renderer": "sequence_conservation", "enabled": enabled,
            "params": {"track_index": track_index}}


FEATURES = {"id": "features", "renderer": "features", "enabled": True, "params": {}}
GC = {"id": "gc_content", "renderer": "dinucleotide_content", "enabled": True, "params": {}}
CASES: tuple[dict[str, Any], ...] = (
    {"id": "no-stack-three-rings", "count": 3, "slots": None},
    {"id": "stack-in-series-order", "count": 3, "slots": [FEATURES, _ring(1), _ring(2), _ring(3), GC]},
    {"id": "stack-reordered", "count": 3, "slots": [FEATURES, _ring(3), _ring(1), _ring(2)]},
    {"id": "track-index-as-text", "count": 3, "slots": [FEATURES, _ring("2"), _ring("1"), _ring("3")]},
    {"id": "disabled-row", "count": 3, "slots": [FEATURES, _ring(2, enabled=False), _ring(3), _ring(1)]},
    {"id": "row-repeated", "count": 3, "slots": [_ring(2), _ring(2, slot_id="ring_2_again"), _ring(1), _ring(3)]},
    {"id": "ring-without-row", "count": 3, "slots": [FEATURES, _ring(3), _ring(1)]},
    {"id": "row-of-no-ring", "count": 2, "slots": [_ring(5), _ring(2), _ring(1)]},
    {"id": "twelve-rings-reversed", "count": 12, "slots": [_ring(index) for index in range(12, 0, -1)]},
)


def _drawn_places(case: dict[str, Any], cfg: GbdrawConfig) -> list[int]:
    tracks = [
        ConservationTrack(source_index=place, track_index=place + 1, track_label=f"ring {place + 1}", track_color=None,
                          reference_side=None, hits=pd.DataFrame())
        for place in range(int(case["count"]))
    ]
    if case["slots"] is None:
        slots = _insert_conservation_slots(
            default_circular_track_slots(),
            _conservation_slots_for_tracks(tracks, cfg=cfg, conservation_ring_width=None, conservation_ring_gap=None),
            axis_index=None,
        )
    else:
        slots = [CircularTrackSlot(id=row["id"], renderer=row["renderer"], enabled=row["enabled"], params=row["params"])
                 for row in case["slots"]]
    return [int(track.track_index) - 1 for track, _color in _drawn_ring_colors("ten", slots, tracks)]


def render() -> str:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    cases = [{**case, "expected": {"drawnPlaces": _drawn_places(case, cfg)}} for case in CASES]
    return json.dumps({"schemaVersion": 1, "cases": cases}, indent=2, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail when the stored vectors are stale")
    args = parser.parse_args()
    text = render()
    if args.check:
        if TARGET.read_text(encoding="utf-8") != text:
            print(f"{TARGET.relative_to(REPO_ROOT)} is stale; run {Path(__file__).name}", file=sys.stderr)
            return 1
        return 0
    TARGET.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
