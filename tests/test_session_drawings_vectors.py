"""The drawing vectors that Python and the Web share (B2-D4, DESIGN §4.2).

``tests/fixtures/sessions/drawings-vectors.json`` pins the drawing ID rule:
``gbdraw.session_drawings.allocate_drawing_id`` here and
``allocateDrawingId`` in ``gbdraw/web/js/services/drawing-view.js``
(``tests/web/drawing-ids.test.mjs``) read the same ``ids`` cases.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.session_drawings import allocate_drawing_id, allocate_session_drawing

VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "sessions" / "drawings-vectors.json").read_text(encoding="utf-8")
)


def test_vectors_declare_their_layout() -> None:
    assert VECTORS["schema"] == 1
    assert set(VECTORS) == {"schema", "ids", "upgrades"}


@pytest.mark.parametrize("case", VECTORS["ids"], ids=lambda case: f"{case['mode']}+{'.'.join(case['taken'])}")
def test_drawing_id_follows_the_shared_vectors(case: dict[str, object]) -> None:
    taken = case["taken"]
    assert isinstance(taken, list)
    assert allocate_drawing_id(str(case["mode"]), taken) == case["id"]


def test_a_second_drawing_of_a_mode_is_numbered() -> None:
    assert allocate_session_drawing("linear", requested_id=None, requested_name=None, taken_ids=["linear"]) == (
        "linear-2",
        "Linear 2",
    )
    assert allocate_session_drawing(
        "linear", requested_id="linear-2", requested_name="Linear 2", taken_ids=["circular", "linear"]
    ) == ("linear-2", "Linear 2")
