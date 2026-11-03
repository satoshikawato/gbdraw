"""The flat Session 50 that tag v0.14.5 writes is read into drawings (D-15-190).

Fixture: ``tests/fixtures/sessions/cli-linear-wrapped.v50.gbdraw-session.json.gz``
(see its ``.provenance.json``), a Linear CLI Session with a wrap length, italic
gene names and a fixed identity scale.
"""

from __future__ import annotations

import copy
import gzip
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session import upgrade_session_document
from gbdraw.session_io import CURRENT_SESSION_VERSION, MODE_SCOPED_SESSION_MAX_VERSION, validate_session

REPO = Path(__file__).parents[1]
FIXTURE = REPO / "tests/fixtures/sessions/cli-linear-wrapped.v50.gbdraw-session.json.gz"


def _released() -> dict[str, Any]:
    return json.loads(gzip.decompress(FIXTURE.read_bytes()))


def _config(drawing: dict[str, Any]) -> dict[str, Any]:
    return drawing["renderRequest"]["diagramOptions"]["config"]


def test_the_released_flat_50_is_a_mode_scoped_version() -> None:
    assert _released()["version"] == MODE_SCOPED_SESSION_MAX_VERSION == 50
    assert CURRENT_SESSION_VERSION == 51


def test_python_reads_the_released_flat_50_into_one_linear_drawing_with_its_values(tmp_path: Path) -> None:
    upgraded = upgrade_session_document(_released(), temporary_directory=tmp_path).document.to_dict()
    assert upgraded["version"] == CURRENT_SESSION_VERSION
    drawings = upgraded["drawings"]
    assert [drawing["mode"] for drawing in drawings] == ["linear"]
    config = _config(drawings[0])
    assert config["canvas"]["linear"]["wrap_length"] == 5000
    assert config["labels"]["linear"]["italic_gene_labels"] is True
    assert config["objects"]["blast_match"]["identity_scale"] == [60.0, 95.0]


def test_replaying_the_released_wrapped_cli_session_draws_wrapped_rows(tmp_path: Path) -> None:
    session_path = tmp_path / "released.gbdraw-session.json"
    session_path.write_bytes(gzip.decompress(FIXTURE.read_bytes()))
    replay = subprocess.run(
        [sys.executable, "-P", "-c", "from gbdraw.cli import main; main()", "linear", "--session", session_path.name,
         "-o", "replayed", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True, env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"},
    )
    assert replay.returncode == 0, replay.stderr
    svg = (tmp_path / "replayed.svg").read_text(encoding="utf-8")
    assert 'data-gbdraw-wrap-row="3"' in svg


def test_a_50_holding_drawings_is_refused_with_a_message() -> None:
    upgraded = upgrade_session_document(_released()).document.to_dict()
    combined = copy.deepcopy(upgraded)
    combined["version"] = 50
    with pytest.raises(ValidationError) as refused:
        validate_session(combined)
    print("old combined 50:", refused.value)
    assert "renderRequest" in str(refused.value)
