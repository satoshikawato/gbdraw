"""Python's side of the shared ticks-anchor vectors (GX-18).

The Web's Radius note for a rendered ticks row converts the band centre in
``runMetadata.trackSlotGeometry`` to the anchor that ``r`` pins
(``tickAnchorRadiusFactor`` in ``gbdraw/web/js/services/track-slot-display.js``,
checked by ``tests/web/track-slot-display.test.mjs``). The stored vectors must
equal what the resolver and the geometry serializer give now, so a change to
the default tick length or to the tick side of a label layout fails here until
the vectors and the Web follow it. Rewrite them with
``python tests/test_circular_tick_anchor_vectors.py``.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.diagrams.circular.assemble import _serialize_circular_track_slot_geometry
from gbdraw.diagrams.circular.radial_layout import resolve_circular_radial_layout
from gbdraw.tracks import CircularTrackSlot, ScalarSpec

VECTORS = Path(__file__).parent / "fixtures" / "circular_tick_anchor_vectors.json"
LAYOUTS = ("label_out_tick_in", "label_in_tick_out", "tick_only", "label_only")
TOTAL_LENGTH = 16_000


class _Feature:
    feature_track_id = 0


def _compute_case(radius_px: float, layout: str, side: str, width_px: float | None) -> dict[str, Any]:
    cfg = GbdrawConfig.from_dict(
        modify_config_dict(
            load_config_toml("gbdraw.data", "config.toml"),
            {"labels.circular.scope": "none", "canvas.circular.track_type": "tuckin"},
        )
    )
    record = SimpleNamespace(seq="N" * TOTAL_LENGTH, id="vector")
    canvas = CircularCanvasConfigurator("vector", CircularRenderProfile(cfg), "none", record)
    canvas.radius = radius_px
    ticks = CircularTrackSlot(
        id="ticks",
        renderer="ticks",
        side=side,
        width=None if width_px is None else ScalarSpec(width_px, "px"),
        params={"tick_label_layout": layout},
    )
    features = CircularTrackSlot(
        id="features",
        renderer="features",
        params={"lane_direction": "outside" if side == "inside" else "inside"},
    )
    slots = [features, ticks] if side == "inside" else [ticks, features]
    radial_layout = resolve_circular_radial_layout(
        total_length=TOTAL_LENGTH,
        canvas_config=canvas,
        slots=slots,
        feature_dict={"a": _Feature()},
    )
    payload = _serialize_circular_track_slot_geometry(
        gb_record=record,
        radial_layout=radial_layout,
        layout_slots=slots,
        base_radius_px=radius_px,
    )["records"][0]
    geometry = next(slot for slot in payload["slots"] if slot["slotId"] == "ticks")
    resolved = next(slot for slot in radial_layout.slots if slot.id == "ticks")
    assert resolved.anchor_radius_px is not None
    return {
        "tickLabelLayout": layout,
        "explicitWidthPx": width_px,
        "axisRadiusPx": payload["axisRadiusPx"],
        "geometry": {key: geometry[key] for key in ("side", "widthPx", "radiusFactor")},
        "anchorFactor": float(resolved.anchor_radius_px) / payload["axisRadiusPx"],
    }


def compute_vectors() -> list[dict[str, Any]]:
    # 390 px: the default length is 0.025 R; 200 px: it is the 6 px floor.
    return [
        _compute_case(radius, layout, side, width)
        for radius, layout, side, width in itertools.product(
            (390.0, 200.0), LAYOUTS, ("inside", "outside"), (None, 7.8)
        )
    ]


def test_tick_anchor_vectors_match_python() -> None:
    stored = json.loads(VECTORS.read_text(encoding="utf-8"))["cases"]
    assert stored == compute_vectors(), (
        "Python's ticks anchor changed; rewrite the vectors with "
        "`python tests/test_circular_tick_anchor_vectors.py` and update tickAnchorRadiusFactor"
    )


if __name__ == "__main__":
    VECTORS.write_text(json.dumps({"cases": compute_vectors()}, indent=1) + "\n", encoding="utf-8")
