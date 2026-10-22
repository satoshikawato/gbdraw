"""The Web lift of a schema 1 or 2 render request against Python's decoder.

`liftLegacyRenderRequest` (gbdraw/web/js/services/legacy-render-request.js)
writes a stored schema 1-2 request as a schema-5 request. The shared vector
`tests/fixtures/legacy-render-request-vectors.json` holds every stored request
of the corpus and the arm cases the corpus does not exercise. Python decodes a
stored request and its lift to the same diagram request, unless the case
names the finding of a deliberate difference (PR-D). Resource tables are
stubbed by their names: the vector holds the requests, not their resources.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import fields, replace
from pathlib import Path
from typing import Any

import pandas
import pytest

import gbdraw.session_request_codec as codec
from gbdraw.analysis.protein_colinearity import OrthogroupResult
from gbdraw.tracks import CircularTrackSlot
from gbdraw.tracks.circular import parse_circular_track_slot

ROOT = Path(__file__).parents[1]
VECTORS = ROOT / "tests/fixtures/legacy-render-request-vectors.json"
# Writes each corpus request and case input with its lift (one JS applier of the case edits).
LIFTS = ROOT / "tests/web/helpers/legacy-request-vectors.mjs"


def _resource_ids(value: Any, found: set[str]) -> set[str]:
    if isinstance(value, dict):
        for key, item in value.items():
            if key.lower().endswith("resourceid") and isinstance(item, str):
                found.add(item)
            _resource_ids(item, found)
    elif isinstance(value, list):
        for item in value:
            _resource_ids(item, found)
    return found


def _table_name(table: Any) -> str:
    return str(table["resource"].iloc[0])


def _slot_reading(slot: Any) -> Any:
    """A Circular slot as the renderer reads its gaps, in pixels.

    The schema-5 decoder keeps a string slot as text, which the renderer
    parses with the strict parser. A legacy spacing in pixels (or of 0) fills
    the gaps a slot leaves unset (`radial_layout.py`, `legacy_spacing.resolve`),
    so it reads as those gaps; a factor of the axis radius stays (OV-798).
    """

    if isinstance(slot, str):
        slot = parse_circular_track_slot(slot)
    spacing = getattr(slot, "legacy_spacing", None)
    fill = None
    if spacing is not None and (
        spacing.unit == "px" or spacing.value == 0
        or (slot.inner_gap_px is not None and slot.outer_gap_px is not None)
    ):
        fill, spacing = spacing.value, None
    gaps = {
        name: None if gap is None else float(gap)
        for name, gap in (
            ("inner_gap_px", slot.inner_gap_px if slot.inner_gap_px is not None else fill),
            ("outer_gap_px", slot.outer_gap_px if slot.outer_gap_px is not None else fill),
        )
    }
    if spacing is None and not getattr(slot, "auto_compress", False) and getattr(slot, "preferred_anchor_radius", None) is None:
        return CircularTrackSlot(**{item.name: getattr(slot, item.name) for item in fields(CircularTrackSlot)} | gaps)
    return replace(slot, legacy_spacing=spacing, **gaps) if hasattr(slot, "legacy_spacing") else replace(slot, **gaps)


def _reading(request: Any) -> str:
    """The decoded request in the one form both schemas can mean.

    Schema 1 stored comparison tables by position, which the renderer draws
    between records n and n + 1 (`record_planning.py`: protein tables by
    `index, index + 1`, BLAST files as adjacent pairs); later schemas store the
    indices. Both read as sorted (query, subject, table) triples. Circular
    slots read as the renderer reads their gaps (`_slot_reading`).
    """

    options = request.options
    comparisons = sorted(
        [
            (item.query_record_index, item.subject_record_index, _table_name(item.matches))
            for item in getattr(options, "linear_comparisons", None) or ()
        ]
        + [(index, index + 1, Path(path).stem) for index, path in enumerate(getattr(options, "blast_files", None) or ())]
        + [
            (index, index + 1, _table_name(table))
            for index, table in enumerate(getattr(options, "protein_comparisons", None) or ())
        ]
    )
    changes: dict[str, Any] = {
        name: None for name in ("linear_comparisons", "blast_files", "protein_comparisons") if hasattr(options, name)
    }
    tracks = getattr(options, "tracks", None)
    slots = getattr(tracks, "circular_track_slots", None)
    if slots:
        changes["tracks"] = replace(tracks, circular_track_slots=tuple(_slot_reading(slot) for slot in slots))
    return repr((replace(request, options=replace(options, **changes)), comparisons))


@pytest.mark.browser
def test_each_vector_request_decodes_as_its_lift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    lifts = json.loads(
        subprocess.run([node, str(LIFTS), str(VECTORS)], check=True, capture_output=True, text=True, cwd=ROOT).stdout
    )
    table = lambda file_path, **_: pandas.DataFrame({"resource": [Path(file_path).stem]})  # noqa: E731
    monkeypatch.setattr(codec, "_read_canonical_table", table)
    monkeypatch.setattr(
        codec, "read_comparison_table",
        lambda handle, *, label: pandas.DataFrame({"resource": [label.split("'")[1]]}),
    )
    monkeypatch.setattr(
        codec, "_read_typed_json_resource",
        lambda resource_id, **kw: (
            OrthogroupResult(orthogroups={resource_id: []}, member_by_protein_id={})
            if kw["value_kind"] == "orthogroupResult"
            else (kw["value_kind"], resource_id)
        ),
    )

    def decoded(request: Any) -> str:
        resources = {}
        for resource_id in _resource_ids(request, set()):
            path = tmp_path / f"{resource_id}.txt"
            path.write_text("")
            resources[resource_id] = str(path)
        try:
            return _reading(codec.decode_canonical_request(
                request, resource_paths=resources, output_directory=tmp_path / "output"
            ))
        except codec.CanonicalRequestCodecError as exc:
            return f"refused: {exc}"

    vectors = json.loads(VECTORS.read_text(encoding="utf-8"))
    assert [lift["id"] for lift in lifts["requests"]] == [request["id"] for request in vectors["requests"]]
    for lift in lifts["requests"]:
        assert (lift["value"]["schema"], lift["repairs"]) == (5, []), lift["id"]
        assert decoded(lift["value"]) == decoded(lift["input"]), lift["id"]
    assert [lift["id"] for lift in lifts["cases"]] == [case["id"] for case in vectors["cases"]]
    for case, lift in zip(vectors["cases"], lifts["cases"]):
        if "refused" in case["expect"]:
            # A form no writer stored: the lift refuses it, and so does Python (D-OD13 F3, D-OD16 F2).
            assert lift.get("refused") == case["expect"]["refused"], case["id"]
            assert case["python"]["refused"] and case["python"]["reason"], case["id"]
            assert decoded(lift["input"]).startswith("refused"), case["id"]
            continue
        stored, lifted = decoded(lift["input"]), decoded(lift["value"])
        if case["python"].get("refused"):
            # A form both runtimes refuse: the lift leaves it to the reader.
            assert case["python"]["reason"], case["id"]
            assert stored.startswith("refused") and lifted.startswith("refused"), case["id"]
            continue
        assert not lifted.startswith("refused"), (case["id"], lifted)
        if case["python"]["equal"]:
            assert lifted == stored, case["id"]
        else:
            # A deliberate difference names its finding; it must still differ.
            assert case["python"]["finding"] and case["python"]["reason"], case["id"]
            assert lifted != stored, case["id"]
