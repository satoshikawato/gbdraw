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
import re
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
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
    RenderOutputRequest,
    SessionDrawingSpec,
    build_session_document,
    load_session_document,
    materialize_session,
    save_session_document,
)
from gbdraw.circular import circular_main
from gbdraw.cli import main as cli_main
from gbdraw.cli_utils.session import render_main
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.session import session_drawing_artifacts
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    DRAWINGS_SESSION_MIN_VERSION,
    write_session_json,
)
from gbdraw.web_support.request_render import render_canonical_web_request
from tests.utils.session_drawing_view import drawing_view, without_drawings_move
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
    assert resaved.to_dict()["ui"]["activeDrawingId"] == source["ui"]["mode"]
    for drawing_id, name in (("circular", "fig_circular"), ("linear", "fig_linear")):
        assert [result["name"] for result in session_drawing_artifacts(resaved, drawing_id).fields["results"]] == [name]

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
            drawing_view(linear_document.to_dict())["renderRequest"],
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
                state={"results": web["results"], "editorState": {"featureCatalog": web["featureCatalog"]}},
            ),
        ]
    )
    path = tmp_path / "batch.gbdraw-session.json"
    write_session_json(path, document.to_dict())

    render_main(["--session", str(path), "-o", str(tmp_path / "out" / "fig")])

    assert _svgs(tmp_path / "out") == ["fig_circular_1.svg", "fig_circular_2.svg", "fig_linear.svg"]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_saved_session_keeps_to_slot_and_draws_the_same_span(tmp_path: Path, mode: str) -> None:
    source = tmp_path / "record.gb"
    SeqIO.write(
        SeqRecord(Seq("ATGCGCAT" * 300), id="record", annotations={"molecule_type": "DNA"}),
        source,
        "genbank",
    )
    table = tmp_path / "regions.tsv"
    table.write_text("set_id\tid\tmark\tstart\tend\nregions\tspan\thighlight\t200\t900\n")
    span = (
        "regions:annotations@set_id=regions,side=overlay,anchor_slot=features,"
        "to_slot=gc_content,layer=underlay,z=-1,show_labels=false"
    )
    if mode == "circular":
        main = circular_main
        slots = ["features:features", "gc_content:dinucleotide_content", span]
        slot_args = [arg for slot in slots for arg in ("--circular_track_slot", slot)]
    else:
        main = linear_main
        slots = ["features:features@side=overlay", "gc_content:dinucleotide_content@side=below", span]
        slot_args = ["--linear_track_axis_index", "0"]
        slot_args += [arg for slot in slots for arg in ("--linear_track_slot", slot)]
    session = tmp_path / "span.gbdraw-session.json"
    main(
        ["--gbk", str(source), "--annotation_table", str(table), *slot_args,
         "-o", str(tmp_path / "first"), "--session_output", str(session)]
    )

    def span_params(node: object) -> list[dict[str, Any]]:
        if isinstance(node, dict):
            own = [node["params"]] if node.get("id") == "regions" and "params" in node else []
            return own + [found for value in node.values() for found in span_params(value)]
        if isinstance(node, list):
            return [found for value in node for found in span_params(value)]
        return []

    saved = span_params(json.loads(session.read_text(encoding="utf-8")))
    assert saved and all(params["to_slot"] == "gc_content" for params in saved)

    render_main(["--session", str(session), "-o", str(tmp_path / "again")])

    first = (tmp_path / "first.svg").read_text(encoding="utf-8")
    assert 'data-gbdraw-annotation-id="span"' in first
    assert (tmp_path / "again.svg").read_text(encoding="utf-8") == first


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


@pytest.mark.parametrize(
    "source",
    ["two-mode-session-46", "two-mode-linear-shown.v48.gbdraw-session.json.gz"],
)
def test_a_cli_resave_of_a_two_mode_project_writes_its_drawings(
    source: str, two_drawings: Path, tmp_path: Path
) -> None:
    path = two_drawings if source == "two-mode-session-46" else SESSION_FIXTURES / source
    render_main(["--session", str(path), "-o", str(tmp_path / "fig"), "--save_session"])

    written = json.loads((tmp_path / "fig.gbdraw-session.json").read_text(encoding="utf-8"))
    assert written["version"] == DRAWINGS_SESSION_MIN_VERSION == CURRENT_SESSION_VERSION
    assert "modes" not in written and "otherModeResult" not in written
    resaved = load_session_document(tmp_path / "fig.gbdraw-session.json")
    assert sorted(drawing.mode for drawing in resaved.drawings) == ["circular", "linear"]


def test_a_fresh_cli_save_and_a_python_api_save_write_drawings(tmp_path: Path) -> None:
    source = tmp_path / "one.gb"
    SeqIO.write(
        [SeqRecord(Seq("ATGCGCAT" * 40), id="one", annotations={"molecule_type": "DNA"})],
        source,
        "genbank",
    )
    circular_main(["--gbk", str(source), "-o", str(tmp_path / "cli"), "--save_session"])
    request = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        output=RenderOutputRequest(output_prefix="api"),
    )
    save_session_document(tmp_path / "api.gbdraw-session.json", request)

    for name, mode in (("cli", "circular"), ("api", "linear")):
        written = json.loads((tmp_path / f"{name}.gbdraw-session.json").read_text(encoding="utf-8"))
        assert written["version"] == DRAWINGS_SESSION_MIN_VERSION == CURRENT_SESSION_VERSION
        assert [drawing["id"] for drawing in written["drawings"]] == [mode]
        assert "renderRequest" not in written and "modes" not in written
        document = load_session_document(tmp_path / f"{name}.gbdraw-session.json")
        assert [drawing.mode for drawing in document.drawings] == [mode]


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
    resaved = drawing_view(load_session_document(tmp_path / "crop.gbdraw-session.json").to_dict(), mode="linear")
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


def test_a_schema_2_linear_record_without_selector_draws_every_record(tmp_path: Path) -> None:
    # OV-149: feature-edits-linear-crop-rc.v33 (request schema 2) with a
    # two-record GenBank file as its second card. As Web Load promotes it, the
    # card without selector or region draws both records, keyed <recordKey>:<n>.
    session = json.loads(
        gzip.decompress((SESSION_FIXTURES / "feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz").read_bytes())
    )
    two_records = (REPO_ROOT / "tests" / "fixtures" / "web_batch_two_records.gb").read_bytes()
    session["resources"]["record-2-genbank"].update(
        data=base64.b64encode(two_records).decode("ascii"), size=len(two_records)
    )
    path = tmp_path / "two-records.v33.json"
    path.write_text(json.dumps(session), encoding="utf-8")

    linear_main(["--session", str(path), "-o", str(tmp_path / "two"), "-f", "svg"])

    card = session["renderRequest"]["records"][1]["recordKey"]
    keys = re.findall(r'data-record-key="([^"]+)"', (tmp_path / "two.svg").read_text(encoding="utf-8"))
    assert sorted(keys) == sorted([session["renderRequest"]["records"][0]["recordKey"], f"{card}:1", f"{card}:2"])


def _cli(monkeypatch: pytest.MonkeyPatch, *args: str) -> None:
    """Run `gbdraw <args>` through the console entry point."""

    monkeypatch.setattr(sys, "argv", ["gbdraw", *args])
    cli_main()


def _fresh_protein_session(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: str = "similarity_groups",
    *extra: str,
) -> Path:
    """A Session saved by a fresh `gbdraw linear` protein render (fake LOSATP)."""

    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from tests.test_render_stage_protein_analysis import _fake_losatp, _request

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    request, _identities = _request(mode, tmp_path)
    sidecar = tmp_path / "fresh.gbdraw-session.json"
    _cli(
        monkeypatch,
        "linear", "--gbk", *(str(record.source.path) for record in request.records), "--losat", "losatp",
        "--losatp_mode", mode, *extra, "-o", str(tmp_path / "fresh"), "-f", "svg",
        "--session_output", str(sidecar),
    )
    return sidecar


def _protein_comparison_kinds(session: dict[str, Any]) -> list[tuple[str, str | None]]:
    return [
        (comparison["kind"], comparison.get("mode"))
        for comparison in session["renderRequest"]["comparisons"]
    ]


# The Web Session writer's form after a protein Generate (session-request.js):
# the pipeline marker reuses the stored result; Similarity groups also stores
# each drawn pair's table, Collinear only its typed result.
_STORED_PROTEIN_FORMS = {
    "similarity_groups": [
        ("precomputedProteinComparison", None),
        ("precomputedProteinComparison", None),
        ("orthogroupResult", None),
        ("generatedProteinComparison", "none"),
    ],
    "collinear": [
        ("collinearityResult", None),
        ("generatedProteinComparison", "none"),
    ],
}


@pytest.mark.parametrize("mode", sorted(_STORED_PROTEIN_FORMS))
def test_a_cli_protein_session_stores_the_results_the_web_stores(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    """OV-1068: a CLI Session stores the resolved results, so its replay searches nothing."""

    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from tests.test_render_stage_protein_analysis import _no_search

    sidecar = _fresh_protein_session(tmp_path, monkeypatch, mode)
    session = drawing_view(json.loads(sidecar.read_text(encoding="utf-8")), mode="linear")
    assert _protein_comparison_kinds(session) == _STORED_PROTEIN_FORMS[mode]
    assert session["losatCache"]["entries"]
    assert session["losatDerivedCache"] == {"entries": []}

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _no_search)
    resaved = tmp_path / "replayed.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(sidecar), "-o", str(tmp_path / "replayed"),
         "-f", "svg", "--session_output", str(resaved))
    assert (tmp_path / "replayed.svg").read_text(encoding="utf-8") == (
        tmp_path / "fresh.svg"
    ).read_text(encoding="utf-8")
    replayed = drawing_view(json.loads(resaved.read_text(encoding="utf-8")), mode="linear")
    assert _protein_comparison_kinds(replayed) == _STORED_PROTEIN_FORMS[mode]


def test_a_cli_protein_session_stores_the_tables_a_gallery_rebuild_stores(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1068: with a reversed record the CLI stores the tables the Gallery tool stores."""

    from tools.refresh_gallery_sessions import (
        _canonical_table_rows,
        _reopen_stored_protein_comparisons,
        _resolved_protein_comparisons,
        _resource_bytes,
    )

    # The Gallery tool reads the mode-scoped layout (its Sessions stay Session 48),
    # so the CLI save is read before its move to drawings.
    with without_drawings_move():
        sidecar = _fresh_protein_session(
            tmp_path, monkeypatch, "similarity_groups",
            "--reverse_complement", "0", "--reverse_complement", "1", "--reverse_complement", "0",
        )
    session = json.loads(sidecar.read_text(encoding="utf-8"))
    assert _protein_comparison_kinds(session) == _STORED_PROTEIN_FORMS["similarity_groups"]
    rebuild = json.loads(json.dumps(session))
    assert _reopen_stored_protein_comparisons(rebuild, "orthogroup") == 1
    resolved = _resolved_protein_comparisons(rebuild)
    assert resolved is not None
    rebuilt, descriptors = resolved

    def stored(source: dict[str, Any], comparisons: list[dict[str, Any]]) -> list[object]:
        return [
            _canonical_table_rows(_resource_bytes(source, comparison["resourceId"]))
            if comparison["kind"] == "precomputedProteinComparison"
            else _resource_bytes(source, comparison["resourceId"])
            for comparison in comparisons
            if comparison["kind"] != "generatedProteinComparison"
        ]

    assert stored(session, session["renderRequest"]["comparisons"]) == stored(
        {"resources": descriptors}, rebuilt
    )


def test_a_fresh_cli_render_builds_no_derived_protein_entries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1050 (D-83): the CLI enters the switch the Web wrapper enters."""

    import gbdraw.api.request_render as request_render_module

    def not_built(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("the CLI built derived entries nobody reads")

    monkeypatch.setattr(request_render_module, "derived_losatp_identity", not_built)
    session = drawing_view(json.loads(_fresh_protein_session(tmp_path, monkeypatch).read_text(encoding="utf-8")), mode="linear")
    assert session["losatCache"]["entries"]
    assert session["losatDerivedCache"] == {"entries": []}


def _session_with_derived_entries(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, list[dict[str, Any]]]:
    """A fresh CLI protein Session whose derived cache holds the Python API's entries."""

    from gbdraw.session import render_session_drawings

    sidecar = _fresh_protein_session(tmp_path, monkeypatch)
    # The Python API, outside the CLI's switch, builds the entries of that Session.
    monkeypatch.chdir(tmp_path)
    with materialize_session(load_session_document(sidecar), output_directory=tmp_path / "api") as materialized:
        (api,) = render_session_drawings(materialized, output_prefix="api", overwrite=True).values()
    built = [dict(entry) for entry in api.losat_derived_cache_entries]
    pairs = [pair for entry in built for pair in entry["payload"]["pairs"]]
    assert pairs and all(isinstance(pair["canonicalTsv"], str) for pair in pairs)
    session = json.loads(sidecar.read_text(encoding="utf-8"))
    (linear,) = (drawing for drawing in session["drawings"] if drawing["mode"] == "linear")
    linear["losatDerivedCache"] = {"entries": built}
    source = tmp_path / "with-derived.gbdraw-session.json"
    source.write_text(json.dumps(session), encoding="utf-8")
    return source, built


def test_a_python_api_replay_that_builds_derived_entries_drops_the_legacy_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1129: outside the CLI's switch the render's own entries supersede the evidence too."""

    from gbdraw.session import render_session_drawings

    sidecar = _fresh_protein_session(tmp_path, monkeypatch)
    session = json.loads(sidecar.read_text(encoding="utf-8"))
    assert drawing_view(session, mode="linear")["losatDerivedCache"] == {"entries": []}
    legacy_derived = {"schema": 1, "kind": "derived-losatp-payload", "key": "legacy-derived", "payload": {"pairs": []}}
    session["legacyArtifacts"] = {"proteinDerivedEvidence": {"schema": 1, "entries": [legacy_derived]}}
    source = tmp_path / "with-evidence.gbdraw-session.json"
    source.write_text(json.dumps(session), encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    with materialize_session(load_session_document(source), output_directory=tmp_path / "api") as materialized:
        (api,) = render_session_drawings(materialized, output_prefix="api", overwrite=True).values()
    assert api.losat_derived_cache_entries
    assert api.migration_report.protein_derived_evidence == ()


def _derived(path: Path) -> list[dict[str, Any]]:
    return drawing_view(load_session_document(path).to_dict(), mode="linear")["losatDerivedCache"]["entries"]


def test_a_cli_replay_keeps_the_derived_entries_of_its_session(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1050 (D-83): a replay builds none and writes the Session's own, unchanged."""

    source, built = _session_with_derived_entries(tmp_path, monkeypatch)
    assert _derived(source) == built
    replayed = tmp_path / "replayed.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(source), "-o", str(tmp_path / "replayed"),
         "-f", "svg", "--session_output", str(replayed))
    assert _derived(replayed) == built


_PROTEIN_ARTIFACTS = ("losatCache", "losatDerivedCache", "proteinIdentityManifest")


def _circular_batch_session_with_protein_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, raw_entries: bool, batch: bool = True
) -> Path:
    """A Circular batch Session (``batch=False``: one Circular diagram) over the records of a protein Session, holding its artifacts.

    A Session keeps its LOSATP artifacts when the drawing changes mode; the
    raw entries, the Python API's derived entries and the manifest are real.
    """

    linear_source, _built = _session_with_derived_entries(tmp_path, monkeypatch)
    linear = drawing_view(json.loads(linear_source.read_text(encoding="utf-8")), mode="linear")
    records = tuple(
        RecordInput(source=GenBankInputSource(tmp_path / f"rec{index}.gbk"))
        for index in range(3)
    )
    request: CircularDiagramRequest | CircularBatchRequest = (
        CircularBatchRequest(records=records, output_policy=CircularBatchOutputPolicy(output_prefix="batch"))
        if batch
        else CircularDiagramRequest(records=records[:1])
    )
    session = build_session_document(request).to_dict()
    (circular,) = session["drawings"]
    for key in _PROTEIN_ARTIFACTS:
        circular[key] = linear[key]
    if not raw_entries:
        circular["losatCache"] = {"entries": []}
    source = tmp_path / "batch-protein.gbdraw-session.json"
    write_session_json(source, session)
    return source


def _protein_artifacts(path: Path) -> dict[str, str]:
    session = drawing_view(load_session_document(path).to_dict(), mode="circular")
    return {key: json.dumps(session[key], sort_keys=True) for key in _PROTEIN_ARTIFACTS}


@pytest.mark.parametrize("raw_entries", [True, False], ids=["raw-and-derived", "derived-only"])
def test_a_circular_batch_replay_keeps_the_protein_artifacts_of_its_session(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, raw_entries: bool
) -> None:
    """OV-1064: the batch re-save keeps the manifest, so its entries still resolve."""

    source = _circular_batch_session_with_protein_artifacts(tmp_path, monkeypatch, raw_entries=raw_entries)
    saved = _protein_artifacts(source)
    assert json.loads(saved["losatDerivedCache"])["entries"]
    assert bool(json.loads(saved["losatCache"])["entries"]) is raw_entries
    replayed = tmp_path / "replayed.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(source), "-o", str(tmp_path / "out" / "replayed"),
         "-f", "svg", "--session_output", str(replayed))
    assert _svgs(tmp_path / "out") == [f"replayed_{index}.svg" for index in (1, 2, 3)]
    assert _protein_artifacts(replayed) == saved


@pytest.mark.parametrize("batch", [False, True], ids=["one-diagram", "batch"])
def test_a_circular_replay_keeps_the_legacy_protein_evidence_of_its_session(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, batch: bool
) -> None:
    """OV-1118: a re-save keeps ``legacyArtifacts`` for one diagram and a batch alike."""

    source = tmp_path / "records.gb"
    SeqIO.write(
        [
            SeqRecord(Seq("ATGCGCAT" * 50), id=f"record{index}", annotations={"molecule_type": "DNA"})
            for index in (1, 2)
        ],
        source,
        "genbank",
    )
    if batch:
        drawing: CircularDiagramRequest | CircularBatchRequest = CircularBatchRequest(
            records=(RecordInput(source=GenBankInputSource(source), cardinality=RecordCardinality.ALL),),
            output_policy=CircularBatchOutputPolicy(output_prefix="batch"),
        )
    else:
        drawing = CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(source), cardinality=RecordCardinality.FIRST),)
        )
    legacy_raw = {"schema": 2, "kind": "raw-losat", "key": "legacy-key", "text": "p_r_old\tp_r_other\n", "program": "blastp"}
    legacy_derived = {"schema": 1, "kind": "derived-losatp-payload", "key": "legacy-derived", "payload": {"pairs": []}}
    # The writer quarantines both legacy entries into legacyArtifacts.
    saved_path = tmp_path / "legacy.gbdraw-session.json"
    save_session_document(
        saved_path,
        drawings=[
            SessionDrawingSpec(
                drawing,
                state={"losatCache": {"entries": [legacy_raw]}, "losatDerivedCache": {"entries": [legacy_derived]}},
            )
        ],
    )
    saved = load_session_document(saved_path).to_dict()["legacyArtifacts"]
    assert set(saved) == {"proteinRawCandidates", "proteinDerivedEvidence"}

    replayed = tmp_path / "replayed.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(saved_path), "-o", str(tmp_path / "out" / "replayed"),
         "-f", "svg", "--session_output", str(replayed))
    assert _svgs(tmp_path / "out") == (["replayed_1.svg", "replayed_2.svg"] if batch else ["replayed.svg"])
    resaved = load_session_document(replayed).to_dict().get("legacyArtifacts")
    assert json.dumps(resaved, sort_keys=True) == json.dumps(saved, sort_keys=True)


@pytest.mark.parametrize("shape", ["one-diagram", "batch", "linear"])
def test_a_replay_drops_only_the_legacy_derived_evidence_that_current_derived_entries_supersede(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, shape: str
) -> None:
    """OV-1119, OV-1129: every shape drops superseded derived evidence and keeps the raw candidates."""

    source = (
        _session_with_derived_entries(tmp_path, monkeypatch)[0]
        if shape == "linear"
        else _circular_batch_session_with_protein_artifacts(
            tmp_path, monkeypatch, raw_entries=True, batch=shape == "batch"
        )
    )
    session = json.loads(source.read_text(encoding="utf-8"))
    mode = "linear" if shape == "linear" else "circular"
    derived_entries = drawing_view(session, mode=mode)["losatDerivedCache"]["entries"]
    assert derived_entries
    legacy_derived = {"schema": 1, "kind": "derived-losatp-payload", "key": "legacy-derived", "payload": {"pairs": []}}
    legacy_raw = {"schema": 2, "kind": "raw-losat", "key": "legacy-key", "text": "p_r_old\tp_r_other\n", "program": "blastp"}
    # A writer that quarantined both kinds; read them back in the on-disk shape.
    legacy_path = tmp_path / "legacy.gbdraw-session.json"
    legacy_genome = tmp_path / "legacy-genome.gb"
    SeqIO.write(
        [SeqRecord(Seq("ATGCGCAT" * 50), id="legacy1", annotations={"molecule_type": "DNA"})],
        legacy_genome,
        "genbank",
    )
    save_session_document(
        legacy_path,
        drawings=[
            SessionDrawingSpec(
                CircularDiagramRequest(
                    records=(RecordInput(source=GenBankInputSource(legacy_genome), cardinality=RecordCardinality.FIRST),)
                ),
                state={"losatCache": {"entries": [legacy_raw]}, "losatDerivedCache": {"entries": [legacy_derived]}},
            )
        ],
    )
    legacy = load_session_document(legacy_path).to_dict()["legacyArtifacts"]
    session["legacyArtifacts"] = legacy
    write_session_json(source, session)
    replayed = tmp_path / "replayed.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(source), "-o", str(tmp_path / "out" / "replayed"),
         "-f", "svg", "--session_output", str(replayed))
    resaved_session = load_session_document(replayed).to_dict()
    # The evidence goes only because the re-save writes current derived entries.
    assert drawing_view(resaved_session, mode=mode)["losatDerivedCache"]["entries"] == derived_entries
    resaved = resaved_session.get("legacyArtifacts")
    assert resaved is not None and set(resaved) == {"proteinRawCandidates"}


_LEGACY_PROTEIN_SESSION = SESSION_FIXTURES / "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz"


def _legacy_protein_session() -> dict[str, Any]:
    """The BGC Session 33: 34 legacy protein raw entries, 8 legacy derived entries, draft Similarity groups."""

    return json.loads(gzip.decompress(_LEGACY_PROTEIN_SESSION.read_bytes()))


def _resaved_legacy_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, session: dict[str, Any], name: str
) -> tuple[list[str], int]:
    """The candidate states and the evidence count of a CLI re-save of ``session``."""

    source = tmp_path / f"{name}.gbdraw-session.json"
    source.write_text(json.dumps(session), encoding="utf-8")
    resaved = tmp_path / f"{name}-resaved.gbdraw-session.json"
    _cli(monkeypatch, "render", "--session", str(source), "-o", str(tmp_path / name / "out"),
         "-f", "svg", "--session_output", str(resaved))
    legacy = load_session_document(resaved).to_dict().get("legacyArtifacts") or {}
    return (
        [entry["state"] for entry in legacy.get("proteinRawCandidates", {}).get("entries", [])],
        len(legacy.get("proteinDerivedEvidence", {}).get("entries", [])),
    )


def test_a_cli_resave_keeps_the_legacy_candidates_its_search_did_not_promote(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1143: the 25 Similarity groups searches promote 25 candidates; 9 stay pending, as after a Web Generate."""

    assert _resaved_legacy_artifacts(tmp_path, monkeypatch, _legacy_protein_session(), "groups") == (
        ["pending"] * 9, 8
    )


def test_a_cli_resave_tries_only_the_candidates_of_the_searches_it_plans(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1143: Pairwise searches the 4 adjacent edges; no protein search tries none."""

    pairwise = _legacy_protein_session()
    pairwise["config"]["losat"]["blastp"]["mode"] = "pairwise"
    assert _resaved_legacy_artifacts(tmp_path, monkeypatch, pairwise, "pairwise") == (["pending"] * 30, 8)
    without = _legacy_protein_session()
    without["renderRequest"]["comparisons"] = []
    assert _resaved_legacy_artifacts(tmp_path, monkeypatch, without, "none") == (["pending"] * 34, 8)


def test_a_non_losatp_search_tries_no_legacy_candidate() -> None:
    """A LOSATN search is not a protein search: no record pair is attempted, so all candidates stay pending."""

    from types import SimpleNamespace

    from gbdraw.api.session_compat import _legacy_candidate_attempt_pairs

    request = SimpleNamespace(options=SimpleNamespace(losat_search=SimpleNamespace(program="losatn")))
    assert _legacy_candidate_attempt_pairs(request, "orthogroup") == ()  # type: ignore[arg-type]


def test_a_cli_resave_that_promotes_the_last_candidate_drops_the_legacy_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1132: with no candidate left pending the evidence goes, as after a Web Generate."""

    from gbdraw.api.session_compat import adapt_session_request
    from gbdraw.session import session_drawing_artifacts, session_to_request

    document = load_session_document(_LEGACY_PROTEIN_SESSION)
    with materialize_session(document, output_directory=tmp_path) as materialized:
        adapted = adapt_session_request(session_to_request(materialized), session_drawing_artifacts(document))
    unmatched = {candidate["originalEntry"]["key"] for candidate in adapted.migration_report.protein_raw_candidates}
    assert len(unmatched) == 9
    session = _legacy_protein_session()
    session["losatCache"]["entries"] = [
        entry for entry in session["losatCache"]["entries"] if entry["key"] not in unmatched
    ]
    assert _resaved_legacy_artifacts(tmp_path, monkeypatch, session, "promoted") == ([], 0)



# Session 49: several drawings of one mode (DESIGN §8 row B2a, B3 / OV-139).
DRAWINGS_V49 = SESSION_FIXTURES / "drawings-v49"


def _linear_pair() -> Path:
    """Drawings ``linear`` (one LOSATN entry) and ``linear-2``, the same request."""

    return DRAWINGS_V49 / "drawings-linear-pair-losat.v49.gbdraw-session.json.gz"


def _circular_pair(tmp_path: Path) -> Path:
    """Drawings ``circular`` and ``circular-2``, each with the committed render
    of ``drawings-circular-only.v49``, written by the Session 49 writer."""

    from gbdraw.session_drawings import write_session_drawings

    only = load_session_document(DRAWINGS_V49 / "drawings-circular-only.v49.gbdraw-session.json.gz")
    view = dict(session_drawing_artifacts(only, "circular").fields)
    path = tmp_path / "circular-pair.v49.gbdraw-session.json"
    write_session_json(path, write_session_drawings([view, view], version=49))
    return path


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_a_mode_command_renders_one_of_several_drawings_of_its_mode(
    mode: str,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from tests.utils.svg_compare import compare_svgs

    session = _linear_pair() if mode == "linear" else _circular_pair(tmp_path)
    main = linear_main if mode == "linear" else circular_main
    name = "Linear" if mode == "linear" else "Circular"

    # Without --drawing the command names every drawing of its mode.
    with pytest.raises(
        ValidationError,
        match=rf"2 {mode} drawings; select one: {mode} \({mode}, '{name}'\), "
        rf"{mode}-2 \({mode}, '{name} 2'\)\.$",
    ):
        main(["--session", str(session), "-o", str(tmp_path / "none"), "-f", "svg"])
    assert _svgs(tmp_path) == []

    render_main(["--session", str(session), "--list_drawings"])
    assert capsys.readouterr().out.splitlines() == [
        f"{mode}\t{mode}\t{name}\tyes",
        f"{mode}-2\t{mode}\t{name} 2\tyes",
    ]

    main(["--session", str(session), "--drawing", f"{mode}-2", "-o", str(tmp_path / "second"), "-f", "svg"])
    assert _svgs(tmp_path) == ["second.svg"]
    # The CLI draws the drawing as the Python API does.
    with materialize_session(load_session_document(session), output_directory=tmp_path / "python") as materialized:
        from gbdraw.api import render_session

        python = render_session(materialized, drawing=f"{mode}-2")
    python_svg = next(path for path in python.output_paths if path.suffix == ".svg")
    result = compare_svgs(tmp_path / "second.svg", python_svg)
    assert result.equal, result.message


def _losat_entries(document: Any) -> dict[str, list[str]]:
    return {
        drawing["id"]: [entry["key"] for entry in drawing.get("losatCache", {}).get("entries", [])]
        for drawing in document.to_dict()["drawings"]
    }


def test_a_resave_keeps_the_losat_entries_of_every_drawing(tmp_path: Path) -> None:
    """B3 / OV-139: each drawing keeps the LOSAT entries it carries, rendered or not."""

    data = json.loads(gzip.decompress(_linear_pair().read_bytes()))
    linear, second = data["drawings"]
    (entry,) = linear["losatCache"]["entries"]
    # linear-2 holds an entry of its own that linear does not hold.
    own = {**entry, "key": "f" * len(entry["key"]), "text": entry["text"].replace("1847", "1846")}
    second["losatCache"] = {"entries": [own]}
    session = tmp_path / "pair.v49.gbdraw-session.json"
    write_session_json(session, data)
    expected = {"linear": [entry["key"]], "linear-2": [own["key"]]}
    assert _losat_entries(load_session_document(session)) == expected

    for selected in ("linear", "linear-2"):
        output = tmp_path / f"{selected}.gbdraw-session.json"
        linear_main(
            ["--session", str(session), "--drawing", selected, "-o", str(tmp_path / selected),
             "-f", "svg", "--session_output", str(output)]
        )
        resaved = load_session_document(output)
        assert resaved.version == 49
        assert _losat_entries(resaved) == expected, selected

    render_main(["--session", str(session), "-o", str(tmp_path / "all"), "-f", "svg", "--save_session"])
    assert _losat_entries(load_session_document(tmp_path / "all.gbdraw-session.json")) == expected
