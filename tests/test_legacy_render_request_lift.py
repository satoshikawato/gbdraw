"""The Web lift of a schema 1 or 2 render request against Python's decoder.

`liftLegacyRenderRequest` (gbdraw/web/js/services/legacy-render-request.js)
writes each stored schema 1-2 request of the corpus as a schema-5 request.
Python decodes the stored request and its lift to the same diagram request
(PR-D). Resource tables are stubbed by their names: the corpus holds the
requests, not their resources.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from dataclasses import replace
from pathlib import Path
from typing import Any

import pandas
import pytest

import gbdraw.session_request_codec as codec
from gbdraw.analysis.protein_colinearity import OrthogroupResult
from gbdraw.tracks.circular import parse_circular_track_slot

ROOT = Path(__file__).parents[1]
CORPUS = ROOT / "tests/fixtures/sessions/legacy-render-requests.vectors.json.gz"
LIFT_CORPUS = """
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
const { liftLegacyRenderRequest } = await import('./gbdraw/web/js/services/legacy-render-request.js');
const corpus = JSON.parse(gunzipSync(readFileSync(process.argv[1])));
process.stdout.write(JSON.stringify(corpus.requests.map(({ id, renderRequest }) => {
  const { value, repairs } = liftLegacyRenderRequest(renderRequest);
  return { id, value, repairs };
})));
"""


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


def _reading(request: Any) -> Any:
    """The decoded request in the one form both schemas can mean.

    A schema-1 protein table is drawn between records n and n + 1, as the
    explicit indices of a later schema say; a string Circular slot is parsed
    when the diagram is drawn, as the strict parser reads it.
    """
    options = request.options
    changes: dict[str, Any] = {}
    tables = getattr(options, "protein_comparisons", None)
    if tables is not None and getattr(options, "linear_comparisons", None) is None:
        changes["protein_comparisons"] = None
        changes["linear_comparisons"] = tuple(
            codec.LinearComparison(index, index + 1, table) for index, table in enumerate(tables)
        )
    tracks = getattr(options, "tracks", None)
    slots = getattr(tracks, "circular_track_slots", None)
    if slots:
        changes["tracks"] = replace(
            tracks,
            circular_track_slots=tuple(
                parse_circular_track_slot(slot) if isinstance(slot, str) else slot for slot in slots
            ),
        )
    return replace(request, options=replace(options, **changes))


@pytest.mark.browser
def test_each_corpus_request_decodes_as_its_lift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")
    lifts = json.loads(
        subprocess.run(
            [node, "--input-type=module", "-e", LIFT_CORPUS, str(CORPUS)],
            check=True, capture_output=True, text=True, cwd=ROOT,
        ).stdout
    )
    monkeypatch.setattr(
        codec, "_read_canonical_table",
        lambda file_path, **_: pandas.DataFrame({"resource": [Path(file_path).name]}),
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
        return repr(_reading(codec.decode_canonical_request(
            request, resource_paths=resources, output_directory=tmp_path / "output"
        )))

    corpus = json.loads(gzip.decompress(CORPUS.read_bytes()))
    assert [lift["id"] for lift in lifts] == [request["id"] for request in corpus["requests"]]
    for stored, lift in zip(corpus["requests"], lifts):
        assert (lift["value"]["schema"], lift["repairs"]) == (5, []), stored["id"]
        assert decoded(lift["value"]) == decoded(stored["renderRequest"]), stored["id"]
