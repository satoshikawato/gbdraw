"""``gbdraw render --session``: render the drawings of a saved Session (PY-D).

The command renders every drawing with a committed render, or the ones named
with ``--drawing``; several drawings write ``<prefix>_<ID>``; every diagram
path and the Session output are checked before the first write; and a re-save
replaces only the rendered drawings of a document validated before any render.
``gbdraw circular|linear --session`` render the Session's drawing of their mode.
"""

from __future__ import annotations

import base64
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
    RenderOutputRequest,
    SessionDrawingSpec,
    build_session_document,
    load_session_document,
    materialize_session,
)
from gbdraw.circular import circular_main
from gbdraw.cli_utils.session import render_main
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.session import session_drawing_artifacts
from gbdraw.session_io import CURRENT_SESSION_VERSION, write_session_json
from gbdraw.web_support.request_render import render_canonical_web_request
from tests.utils.two_mode_session import two_mode_session

REPO_ROOT = Path(__file__).resolve().parents[1]
SESSION_FIXTURES = REPO_ROOT / "tests" / "fixtures" / "sessions"


@pytest.fixture
def two_drawings(tmp_path: Path) -> Path:
    """A Session 46 with a Circular and a Linear drawing (``otherModeResult``)."""

    path = tmp_path / "two.gbdraw-session.json"
    write_session_json(path, two_mode_session())
    return path


def _svgs(directory: Path) -> list[str]:
    return sorted(path.name for path in directory.glob("*.svg"))


def test_list_drawings_prints_one_line_per_drawing(
    two_drawings: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    render_main(["--session", str(two_drawings), "--list_drawings"])

    assert capsys.readouterr().out.splitlines() == [
        "circular\tcircular\tCircular\tyes",
        "linear\tlinear\tLinear\tyes",
    ]
    settings_only = SESSION_FIXTURES / "settings-only.v42.json.gz"
    render_main(["--session", str(settings_only), "--list_drawings"])
    assert capsys.readouterr().out.splitlines() == ["circular\tcircular\tCircular\tno"]


def test_render_writes_each_drawing_under_its_id(two_drawings: Path, tmp_path: Path) -> None:
    render_main(["--session", str(two_drawings), "-o", str(tmp_path / "out" / "fig"), "-f", "svg"])

    assert _svgs(tmp_path / "out") == ["fig_circular.svg", "fig_linear.svg"]


def test_render_without_output_uses_each_saved_prefix(
    two_drawings: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    session = two_mode_session()
    prefixes = (
        session["renderRequest"]["output"]["prefix"],
        session["otherModeResult"]["renderRequest"]["output"]["prefix"],
    )

    render_main(["--session", str(two_drawings), "-f", "svg"])

    assert _svgs(work) == sorted(
        (f"{prefixes[0]}_circular.svg", f"{prefixes[1]}_linear.svg")
    )


def test_one_named_drawing_keeps_its_names(two_drawings: Path, tmp_path: Path) -> None:
    render_main(["--session", str(two_drawings), "--drawing", "Linear", "-o", str(tmp_path / "one")])
    render_main(
        ["--session", str(two_drawings), "--drawing", "linear", "--drawing", "circular",
         "-o", str(tmp_path / "both" / "fig")]
    )

    assert _svgs(tmp_path) == ["one.svg"]
    assert _svgs(tmp_path / "both") == ["fig_circular.svg", "fig_linear.svg"]


def test_every_output_is_checked_before_the_first_write(two_drawings: Path, tmp_path: Path) -> None:
    occupied = tmp_path / "fig_linear.svg"
    occupied.write_text("keep", encoding="utf-8")

    with pytest.raises(ValidationError, match="already exist"):
        render_main(["--session", str(two_drawings), "-o", str(tmp_path / "fig")])
    assert _svgs(tmp_path) == ["fig_linear.svg"]
    assert occupied.read_text(encoding="utf-8") == "keep"

    sidecar = tmp_path / "fig.gbdraw-session.json"
    sidecar.write_text("keep", encoding="utf-8")
    occupied.unlink()
    with pytest.raises(ValidationError, match="Session output already exists"):
        render_main(["--session", str(two_drawings), "-o", str(tmp_path / "fig"), "--save_session"])
    assert _svgs(tmp_path) == []
    assert sidecar.read_text(encoding="utf-8") == "keep"


def test_save_session_replaces_the_rendered_drawings(two_drawings: Path, tmp_path: Path) -> None:
    source = two_mode_session()

    with pytest.raises(ValidationError, match="several drawings needs -o"):
        render_main(["--session", str(two_drawings), "--save_session"])
    with pytest.raises(SystemExit):
        render_main(
            ["--session", str(two_drawings), "--save_session", "--session_output", "x.json"]
        )

    render_main(["--session", str(two_drawings), "-o", str(tmp_path / "fig"), "--save_session"])
    resaved = load_session_document(tmp_path / "fig.gbdraw-session.json")
    assert [drawing.id for drawing in resaved.drawings] == ["circular", "linear"]
    assert resaved.to_dict()["ui"]["mode"] == source["ui"]["mode"]
    assert [result["name"] for result in resaved.to_dict()["results"]] == ["fig_circular"]
    assert [result["name"] for result in resaved.to_dict()["otherModeResult"]["results"]] == ["fig_linear"]

    # Only the named drawing is replaced; the other keeps its saved set.
    output = tmp_path / "linear-only.gbdraw-session.json.gz"
    render_main(
        ["--session", str(two_drawings), "--drawing", "linear", "-o", str(tmp_path / "lin"),
         "--session_output", str(output)]
    )
    partial = load_session_document(output)
    circular = session_drawing_artifacts(partial, "circular").fields
    assert circular["renderRequest"] == source["renderRequest"]
    assert circular["results"] == source["results"]
    linear = session_drawing_artifacts(partial, "linear").fields
    assert [result["name"] for result in linear["results"]] == ["lin"]


def test_a_circular_batch_drawing_still_numbers_its_diagrams(tmp_path: Path) -> None:
    source = tmp_path / "two-records.gb"
    SeqIO.write(
        [
            SeqRecord(Seq("ATGCGCAT" * 40), id=record_id, annotations={"molecule_type": "DNA"})
            for record_id in ("first", "second")
        ],
        source,
        "genbank",
    )
    record = RecordInput(source=GenBankInputSource(source), cardinality=RecordCardinality.ALL)
    linear = LinearDiagramRequest(records=(record,), output=RenderOutputRequest(output_prefix="lin"))
    linear_document = build_session_document(linear)
    with materialize_session(linear_document, output_directory=tmp_path / "web") as materialized:
        web = render_canonical_web_request(
            linear_document.to_dict()["renderRequest"],
            resource_paths=materialized.resource_paths,
            output_directory=tmp_path / "web-output",
        )
    document = build_session_document(
        drawings=[
            CircularBatchRequest(
                records=(record,),
                output_policy=CircularBatchOutputPolicy(output_prefix="batch"),
            ),
            SessionDrawingSpec(
                linear,
                state={"results": web["results"], "editorState": {"featureCatalog": web["metadata"]["featureCatalog"]}},
            ),
        ]
    )
    path = tmp_path / "batch.gbdraw-session.json"
    write_session_json(path, document.to_dict())

    render_main(["--session", str(path), "-o", str(tmp_path / "out" / "fig")])

    assert _svgs(tmp_path / "out") == ["fig_circular_1.svg", "fig_circular_2.svg", "fig_linear.svg"]


def test_a_session_without_a_committed_render_cannot_render(tmp_path: Path) -> None:
    settings_only = SESSION_FIXTURES / "settings-only.v42.json.gz"

    with pytest.raises(ValidationError, match="Settings-only Session"):
        render_main(["--session", str(settings_only), "-o", str(tmp_path / "out")])
    assert _svgs(tmp_path) == []


def test_render_replays_a_session_27_30_with_its_mode(tmp_path: Path) -> None:
    source_path = tmp_path / "legacy.gb"
    SeqIO.write(
        SeqRecord(Seq("ATGCGCAT"), id="legacy", annotations={"molecule_type": "DNA"}),
        source_path,
        "genbank",
    )
    data = source_path.read_bytes()
    session = {
        "format": "gbdraw-session",
        "version": 30,
        "createdAt": "2026-07-30T00:00:00Z",
        "config": {"form": {"prefix": "legacy"}, "adv": {}},
        "ui": {"mode": "linear", "lInputType": "gb"},
        "files": {"linearSeqs": [{"gb": {
            "name": source_path.name, "type": "application/octet-stream", "size": len(data),
            "lastModified": 0, "data": base64.b64encode(data).decode("ascii"),
        }}]},
        "cliInvocation": {
            "schema": 1, "mode": "linear",
            "args": ["--gbk", source_path.name, "--format", "svg", "--no-gc", "--no-skew"],
            "renderFormats": ["svg"],
            "fileBindings": [{"argIndex": 1, "slot": "files.linearSeqs[0].gb", "name": source_path.name}],
            "generatedBy": "gbdraw",
        },
    }
    path = tmp_path / "legacy.gbdraw-session.json"
    path.write_text(json.dumps(session), encoding="utf-8")

    with pytest.raises(ValidationError, match="no drawing 'circular'"):
        render_main(["--session", str(path), "--drawing", "circular"])
    render_main(["--session", str(path), "-o", str(tmp_path / "out" / "legacy")])

    assert _svgs(tmp_path / "out") == ["legacy.svg"]


def test_mode_commands_render_their_own_drawing(two_drawings: Path, tmp_path: Path) -> None:
    linear_main(["--session", str(two_drawings), "-o", str(tmp_path / "lin"), "-f", "svg"])
    circular_main(["--session", str(two_drawings), "--drawing", "circular", "-o", str(tmp_path / "circ")])

    assert _svgs(tmp_path) == ["circ.svg", "lin.svg"]
    with pytest.raises(ValidationError, match="is a linear drawing, not circular.*circular .*linear"):
        circular_main(["--session", str(two_drawings), "--drawing", "linear", "-o", str(tmp_path / "x")])
    with pytest.raises(SystemExit):
        linear_main(
            ["--session", str(two_drawings), "--drawing", "linear", "--drawing", "circular",
             "-o", str(tmp_path / "y")]
        )
    with pytest.raises(ValidationError, match="--drawing selects a drawing of a --session file"):
        linear_main(["--gbk", "missing.gb", "--drawing", "linear"])
    assert _svgs(tmp_path) == ["circ.svg", "lin.svg"]


def test_a_session_that_cannot_be_saved_again_writes_no_diagram(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_migration as session_migration

    session = json.loads(
        gzip.decompress((SESSION_FIXTURES / "feature-edits-crop-rc.v44.gbdraw-session.json.gz").read_bytes())
    )
    vector = json.loads((REPO_ROOT / "tests" / "fixtures" / "feature-placement-migration.json").read_text(encoding="utf-8"))
    session["config"]["featurePlacementOverrides"] = {
        key: row
        for key, row in vector["input"].items()
        if row["placement"]["kind"] == "main" or row["placement"]["side"] in {"above", "below"}
    }
    path = tmp_path / "placements.v44.json"
    path.write_text(json.dumps(session), encoding="utf-8")
    split = session_migration.split_draft_into_modes

    def split_without_scopes(draft: Any, **context: Any) -> dict[str, Any]:
        # A split that puts every row in the Circular slice: the re-saved
        # Session would be invalid, so the run fails before any diagram exists.
        result = split(draft, **context)
        result["modes"]["circular"]["config"]["featurePlacementOverrides"] = {
            key: row
            for slice_ in result["modes"].values()
            for key, row in slice_["config"]["featurePlacementOverrides"].items()
        }
        return result

    monkeypatch.setattr(session_migration, "split_draft_into_modes", split_without_scopes)

    with pytest.raises(ValidationError, match="unsupported in circular mode"):
        render_main(["--session", str(path), "-o", str(tmp_path / "out"), "--save_session"])
    assert _svgs(tmp_path) == []
    assert not (tmp_path / "out.gbdraw-session.json").exists()


def test_the_resaved_session_is_current(two_drawings: Path, tmp_path: Path) -> None:
    render_main(
        ["--session", str(SESSION_FIXTURES / "q-frame-main-linear-reverse.v42.gbdraw-session.json.gz"),
         "-o", str(tmp_path / "frame"), "--save_session"]
    )

    assert load_session_document(tmp_path / "frame.gbdraw-session.json").version == CURRENT_SESSION_VERSION


def test_a_resave_prints_the_results_the_upgrade_drops(tmp_path: Path) -> None:
    # Session 33 saved no feature catalog: the upgrade drops its Result and
    # names it, and the re-save holds the Result this render writes.
    source = SESSION_FIXTURES / "feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz"
    environment = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    result = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "render", "--session", str(source),
         "-o", str(tmp_path / "crop"), "-f", "svg", "--save_session"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=240,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (
        f"WARNING: Upgrading Session 33 to {CURRENT_SESSION_VERSION} dropped the linear drawing's "
        "1 Result(s) 'out.svg': Session 33 saved no feature catalog for them. "
        "Render the drawing and save the Session to write new Results."
    ) in result.stdout.splitlines()
    resaved = load_session_document(tmp_path / "crop.gbdraw-session.json").to_dict()
    assert [entry["name"] for entry in resaved["results"]] == ["crop"]


def test_a_resave_reads_session_39_table_rows_as_the_current_writer_writes_them(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # OV-148: the Session 39 Label whitelist row has four cells (a tab typed in
    # the keyword). It is read as the current writer writes it, as Web Load reads
    # it, and the re-save holds the table it read, so the 46 reads it again.
    source = SESSION_FIXTURES / "whitelist-tab-keyword.v39.gbdraw-session.json.gz"
    render_main(["--session", str(source), "-o", str(tmp_path / "fig"), "-f", "svg", "--save_session"])

    assert _svgs(tmp_path) == ["fig.svg"]
    assert "cytochrome c oxidase" in (tmp_path / "fig.svg").read_text(encoding="utf-8")
    assert "line(s) 1 had extra cells, joined into the last column with one space" in caplog.text
    resaved = load_session_document(tmp_path / "fig.gbdraw-session.json").to_dict()
    assert resaved["version"] == CURRENT_SESSION_VERSION
    assert base64.b64decode(resaved["resources"]["label-whitelist-table"]["data"]).decode() == (
        "feature_type\tqualifier\tkeyword\nCDS\tproduct\tcytochrome c oxidase\n"
    )
    caplog.clear()
    render_main(["--session", str(tmp_path / "fig.gbdraw-session.json"), "-o", str(tmp_path / "again"), "-f", "svg"])
    assert (tmp_path / "again.svg").read_bytes() == (tmp_path / "fig.svg").read_bytes()
    assert "had extra cells" not in caplog.text
