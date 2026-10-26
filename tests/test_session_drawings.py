"""The drawings of a Session: selection, render-all, writers and the upgrade.

A drawing is one diagram of a Session (project). ``gbdraw.session_drawings``
is the one owner of where a drawing's parts sit in the document; the API, the
CLI and the writers reach drawings only through it. Session 46 holds at most
one drawing of each mode (the second in ``otherModeResult``), named by its
mode; each drawing's draft is its mode's slice of ``modes``.
"""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path

import pytest

import gbdraw.api.request_render as request_render_module
from gbdraw.api import (
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    SessionConversionError,
    SessionDocument,
    SessionDrawing,
    SessionDrawingSelectionError,
    SessionDrawingSpec,
    SessionFormatError,
    SessionRenderError,
    SessionResourceError,
    SessionUpgrade,
    SessionVersionError,
    build_session_document,
    load_session_document,
    materialize_session,
    render_session,
    render_session_drawings,
    save_session_document,
    session_to_request,
    upgrade_session_document,
)
from gbdraw.exceptions import ValidationError
from gbdraw.session_drawings import (
    _view,
    replace_session_drawings,
    write_session_drawings,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION, MODE_SCOPED_SESSION_MAX_VERSION
from gbdraw.web_support.request_render import render_canonical_web_request
from tests.utils.session_drawing_view import drawing_view

REPO_ROOT = Path(__file__).resolve().parents[1]
GENBANK = REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
SESSION_FIXTURES = REPO_ROOT / "tests" / "fixtures" / "sessions"


def _fixture(name: str) -> dict:
    data = (SESSION_FIXTURES / name).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return json.loads(data)


def _requests(prefix: str = "shared") -> tuple[CircularDiagramRequest, LinearDiagramRequest]:
    record = RecordInput(source=GenBankInputSource(GENBANK))
    output = RenderOutputRequest(output_prefix=prefix)
    return (
        CircularDiagramRequest(records=(record,), output=output),
        LinearDiagramRequest(records=(record,), output=output),
    )


def _linear_result_state(request: LinearDiagramRequest, tmp_path: Path) -> dict:
    """The Web-owned Result set of a rendered Linear drawing."""

    document = build_session_document(request)
    with materialize_session(document, output_directory=tmp_path / "web") as materialized:
        web = render_canonical_web_request(
            drawing_view(document.to_dict(), mode="linear")["renderRequest"],
            resource_paths=materialized.resource_paths,
            output_directory=tmp_path / "web-output",
        )
    return {
        "results": web["results"],
        "editorState": {"featureCatalog": web["featureCatalog"]},
    }


def _two_drawing_session(tmp_path: Path, **kwargs) -> SessionDocument:
    """A Circular and a Linear drawing of one GenBank file."""

    circular, linear = _requests()
    return build_session_document(
        drawings=[
            circular,
            SessionDrawingSpec(linear, state=_linear_result_state(linear, tmp_path)),
        ],
        **kwargs,
    )


@pytest.mark.parametrize("mode", ("circular", "linear"))
def test_a_one_request_session_is_one_drawing_named_by_its_mode(mode: str) -> None:
    circular, linear = _requests()
    document = build_session_document(circular if mode == "circular" else linear)

    drawing = SessionDrawing(mode, mode.capitalize(), mode, True)
    assert document.drawings == (drawing,)
    assert document.drawing() == document.drawing(mode) == drawing
    assert document.drawing(mode.capitalize()) == document.drawing(mode=mode) == drawing
    assert document.active_drawing_id == mode
    assert document.mode == mode
    assert document.has_canonical_request is True
    with pytest.raises(
        SessionDrawingSelectionError,
        match=rf"no drawing 'linear-2'; it has: {mode} \({mode}, '{mode.capitalize()}'\)",
    ):
        document.drawing("linear-2")
    other = "linear" if mode == "circular" else "circular"
    with pytest.raises(SessionDrawingSelectionError, match=f"no {other} drawing"):
        document.drawing(mode=other)


def test_a_settings_only_session_is_one_drawing_without_a_render(tmp_path: Path) -> None:
    document = load_session_document(_fixture("settings-only.v42.json.gz"))

    assert document.drawings == (SessionDrawing("circular", "Circular", "circular", False),)
    assert document.mode == "circular"
    assert document.has_canonical_request is False
    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionConversionError, match="Settings-only Session"):
            session_to_request(materialized)
        with pytest.raises(SessionDrawingSelectionError, match="no drawing with a committed render"):
            render_session_drawings(materialized)
        with pytest.raises(SessionDrawingSelectionError, match="'circular' has no committed render"):
            render_session_drawings(materialized, drawings=["circular"])


def test_a_session_27_30_is_one_drawing_of_its_declared_mode() -> None:
    document = load_session_document(_fixture("cli-linear-protein.v30.gbdraw-session.json.gz"))

    assert document.drawings == (SessionDrawing("linear", "Linear", "linear", False),)
    assert document.active_drawing_id == "linear"


def test_render_session_drawings_names_each_drawing_output(tmp_path: Path) -> None:
    document = _two_drawing_session(tmp_path)
    assert [drawing.id for drawing in document.drawings] == ["circular", "linear"]
    assert list(document.to_dict()["resources"]) == ["record-1-genbank"]

    with materialize_session(document, output_directory=tmp_path / "all") as materialized:
        results = render_session_drawings(materialized)
    assert list(results) == ["circular", "linear"]
    assert [path.name for path in results["circular"].output_paths] == ["shared_circular.svg"]
    assert [path.name for path in results["linear"].output_paths] == ["shared_linear.svg"]

    with materialize_session(document, output_directory=tmp_path / "prefixed") as materialized:
        results = render_session_drawings(
            materialized,
            drawings=["Linear", "circular"],
            output_prefix="figure",
            formats="svg",
        )
    # Results follow the document order, not the selector order.
    assert list(results) == ["circular", "linear"]
    assert results["linear"].output_paths == (tmp_path / "prefixed" / "figure_linear.svg",)

    # One drawing keeps today's names.
    with materialize_session(document, output_directory=tmp_path / "one") as materialized:
        (only,) = render_session_drawings(materialized, drawings=["linear"]).values()
    assert only.output_paths == (tmp_path / "one" / "shared.svg",)


def test_render_session_drawings_checks_every_output_before_writing(tmp_path: Path) -> None:
    document = _two_drawing_session(tmp_path)
    occupied = tmp_path / "out" / "shared_linear.svg"
    occupied.parent.mkdir()
    occupied.write_text("keep this diagram", encoding="utf-8")

    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        with pytest.raises(SessionRenderError, match="already exist"):
            render_session_drawings(materialized)
        assert not (tmp_path / "out" / "shared_circular.svg").exists()
        assert occupied.read_text(encoding="utf-8") == "keep this diagram"
        results = render_session_drawings(materialized, overwrite=True)
    assert results["linear"].output_paths == (occupied,)
    assert occupied.read_text(encoding="utf-8").startswith("<?xml")


def test_render_session_drawings_parses_each_resource_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    document = _two_drawing_session(tmp_path)
    original_loader = request_render_module.load_gbks
    loads: list[tuple[str, ...]] = []

    def counting_loader(paths: list[str]):
        loads.append(tuple(paths))
        return original_loader(paths)

    monkeypatch.setattr(request_render_module, "load_gbks", counting_loader)

    with materialize_session(document, output_directory=tmp_path / "together") as materialized:
        render_session_drawings(materialized)
    assert len(loads) == 1

    loads.clear()
    for drawing in ("circular", "linear"):
        with materialize_session(document, output_directory=tmp_path / drawing) as materialized:
            render_session(materialized, drawing=drawing)
    assert len(loads) == 2


def test_drawings_that_share_a_resource_draw_as_each_drawing_alone(tmp_path: Path) -> None:
    """One call parses a shared record once for every drawing; no stage changes the shared parsed records."""

    document = _two_drawing_session(tmp_path)
    with materialize_session(document, output_directory=tmp_path / "together") as materialized:
        together = render_session_drawings(materialized, formats="svg")
    for drawing, result in together.items():
        (path,) = result.output_paths
        with materialize_session(document, output_directory=tmp_path / drawing) as materialized:
            (alone,) = render_session_drawings(
                materialized, drawings=[drawing], output_prefix=path.stem, formats="svg"
            ).values()
        assert alone.output_paths[0].name == path.name
        assert alone.output_paths[0].read_bytes() == path.read_bytes()


def test_session_renders_keep_their_error_surface(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from gbdraw.api.prepared import PreparedBiologicalInputCache
    from gbdraw.cli_utils.session import render_main
    from gbdraw.render.output_paths import OutputExistsError

    document = _two_drawing_session(tmp_path)
    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        pass
    with pytest.raises(SessionResourceError, match="no longer active"):
        render_session(materialized, drawing="circular")
    with pytest.raises(SessionResourceError, match="no longer active"):
        render_session_drawings(materialized)

    # An error of the body keeps its type: the CLI words this refusal with --overwrite.
    source = tmp_path / "two.gbdraw-session.json"
    source.write_text(json.dumps(document.to_dict()), encoding="utf-8")
    sidecar = tmp_path / "existing.gbdraw-session.json"
    sidecar.write_text("keep", encoding="utf-8")
    with pytest.raises(OutputExistsError) as refused:
        render_main(["--session", str(source), "-o", str(tmp_path / "fig"), "--session_output", str(sidecar)])
    assert type(refused.value) is OutputExistsError
    assert refused.value.paths == (str(sidecar),)
    assert sidecar.read_text(encoding="utf-8") == "keep"

    # Opening the cache fails inside another transaction.
    with materialize_session(document, output_directory=tmp_path / "nested") as materialized:
        with PreparedBiologicalInputCache().transaction(resource_paths={}, diagnostics=None):
            with pytest.raises(SessionRenderError, match="Canonical session rendering failed: .*cannot be nested"):
                render_session(materialized, drawing="circular")

    def failing_commit(self, transaction):
        raise ValidationError("probe")

    monkeypatch.setattr(PreparedBiologicalInputCache, "_commit", failing_commit)
    with materialize_session(document, output_directory=tmp_path / "commit") as materialized:
        with pytest.raises(SessionRenderError, match="Canonical session rendering failed: probe"):
            render_session(materialized, drawing="circular")
        with pytest.raises(SessionRenderError, match="Canonical session rendering failed: probe"):
            render_session_drawings(materialized)


def test_drawing_selection_errors_name_the_drawings(tmp_path: Path) -> None:
    document = _two_drawing_session(tmp_path)

    with pytest.raises(SessionDrawingSelectionError, match=r"2 drawings; select one: circular .*linear"):
        document.drawing()
    with pytest.raises(SessionDrawingSelectionError, match="select one"):
        document.mode
    with materialize_session(document, output_directory=tmp_path / "out") as materialized:
        with pytest.raises(SessionDrawingSelectionError, match="more than once: linear"):
            render_session_drawings(materialized, drawings=["linear", "Linear"])
        with pytest.raises(SessionDrawingSelectionError, match="at least one"):
            render_session_drawings(materialized, drawings=[])
        with pytest.raises(SessionDrawingSelectionError, match="not one string"):
            render_session_drawings(materialized, drawings="linear")
        with pytest.raises(SessionDrawingSelectionError, match="no drawing 'linear-2'"):
            render_session_drawings(materialized, drawings=["linear-2"])
        with pytest.raises(SessionDrawingSelectionError):
            render_session(materialized)


def test_build_session_document_writes_two_drawings(tmp_path: Path) -> None:
    document = _two_drawing_session(tmp_path, active_drawing="Linear")

    data = document.to_dict()
    # The drawings keep the given order; ``active_drawing`` only names the one that opens.
    assert [drawing["id"] for drawing in data["drawings"]] == ["circular", "linear"]
    assert drawing_view(data, "circular")["renderRequest"]["mode"] == "circular"
    assert drawing_view(data, "linear")["renderRequest"]["mode"] == "linear"
    assert data["ui"] == {"activeDrawingId": "linear"}
    assert document.active_drawing_id == "linear"

    path = tmp_path / "two.gbdraw-session.json"
    circular, linear = _requests()
    save_session_document(
        path,
        drawings=[
            circular,
            SessionDrawingSpec(linear, state=_linear_result_state(linear, tmp_path / "save")),
        ],
    )
    assert load_session_document(path).drawings == document.drawings


def test_build_session_document_admits_what_session_49_holds() -> None:
    circular, linear = _requests()
    document = build_session_document(
        drawings=[
            circular,
            circular,
            SessionDrawingSpec(mode="linear"),
            SessionDrawingSpec(linear, state={"cliOptions": {"x": 1}}),
        ],
        active_drawing="Linear 2",
    )

    data = document.to_dict()
    assert [(drawing["id"], drawing["name"]) for drawing in data["drawings"]] == [
        ("circular", "Circular"),
        ("circular-2", "Circular 2"),
        ("linear", "Linear"),
        ("linear-2", "Linear 2"),
    ]
    assert data["drawings"][2]["renderRequest"] is None
    # A drawing field belongs to its drawing; another drawing does not share it.
    assert data["drawings"][3]["cliOptions"] == {"x": 1}
    assert "cliOptions" not in data["drawings"][0]
    assert data["ui"] == {"activeDrawingId": "linear-2"}


def test_build_session_document_refuses_drawings_it_cannot_name() -> None:
    circular, linear = _requests()

    with pytest.raises(SessionFormatError, match="names each drawing by its mode"):
        build_session_document(drawings=[SessionDrawingSpec(circular, id="circular-2")])
    with pytest.raises(SessionFormatError, match="names each drawing by its mode"):
        build_session_document(drawings=[SessionDrawingSpec(circular, name="Overview")])
    with pytest.raises(SessionFormatError, match="does not match"):
        build_session_document(drawings=[SessionDrawingSpec(circular, mode="linear")])
    with pytest.raises(SessionFormatError, match="needs a mode"):
        build_session_document(drawings=[SessionDrawingSpec()])
    with pytest.raises(SessionFormatError, match="either request or drawings"):
        build_session_document(circular, drawings=[linear])
    with pytest.raises(SessionFormatError, match="either request or drawings"):
        build_session_document()
    with pytest.raises(SessionDrawingSelectionError, match="names none"):
        build_session_document(circular, active_drawing="linear")


def test_a_drawing_state_cannot_hold_another_modes_slice() -> None:
    """Session 49 retired the shared ``modes`` key: another mode's slice is
    refused, neither dropped nor turned into a drawing of that mode."""

    circular, _linear = _requests()
    other = {"linear": {"config": {"form": {"align_center": True}}}}

    with pytest.raises(ValidationError, match=r"modes\.linear"):
        build_session_document(drawings=[SessionDrawingSpec(circular, state={"modes": other})])
    with pytest.raises(ValidationError, match=r"modes\.linear"):
        build_session_document(
            drawings=[SessionDrawingSpec(circular, state={"modes": {"circular": {}, **other}})]
        )


@pytest.mark.parametrize(
    ("web_files", "field"),
    [
        ({"bindings": {"schema": 2, "linearSeqs": []}}, r"webFiles\.bindings\.linearSeqs"),
        ({"linearRecordMetadata": []}, r"webFiles\.linearRecordMetadata"),
    ],
)
def test_a_drawing_state_cannot_hold_another_modes_web_files(web_files: dict, field: str) -> None:
    """Another mode's Web file bindings and request metadata are refused by
    name, as its ``modes`` slice is, never dropped silently (D-15-125)."""

    circular, _linear = _requests()

    with pytest.raises(ValidationError, match=field):
        build_session_document(drawings=[SessionDrawingSpec(circular, state={"webFiles": web_files})])


def test_each_drawing_keeps_its_own_losat_artifacts() -> None:
    def view(mode: str, entries: list) -> dict:
        return {
            "renderRequest": {"mode": mode},
            "results": [{"name": mode, "content": "<svg/>"}],
            "losatCache": {"entries": entries},
        }

    data = write_session_drawings(
        [view("circular", [{"key": "ring"}]), view("linear", [{"key": "ring"}, {"key": "protein"}])]
    )
    assert [drawing["losatCache"]["entries"] for drawing in data["drawings"]] == [
        [{"key": "ring"}],
        [{"key": "ring"}, {"key": "protein"}],
    ]


def test_each_drawing_gives_its_own_mode_slice(tmp_path: Path) -> None:
    circular, linear = _requests()
    circular_slice = {"config": {"adv": {"axis_stroke_width": 2}}}
    linear_slice = {"config": {"form": {"align_center": True}}}
    document = build_session_document(
        drawings=[
            SessionDrawingSpec(circular, state={"modes": {"circular": circular_slice}}),
            SessionDrawingSpec(
                linear,
                state={**_linear_result_state(linear, tmp_path), "modes": {"linear": linear_slice}},
            ),
        ]
    )
    data = document.to_dict()
    assert drawing_view(data, "circular")["modes"] == {"circular": circular_slice}
    assert drawing_view(data, "linear")["modes"] == {"linear": linear_slice}
    # A Session 46-48 moves to drawings before it is saved again.
    with pytest.raises(ValidationError, match="Only a Session 49"):
        replace_session_drawings(
            {**data, "version": MODE_SCOPED_SESSION_MAX_VERSION},
            {"circular": drawing_view(data, "circular")},
        )


def test_upgrade_returns_a_current_document_unchanged() -> None:
    circular, _linear = _requests()
    document = build_session_document(circular)

    assert upgrade_session_document(document) == SessionUpgrade(document)
    assert upgrade_session_document(document).document is document


def test_upgrade_refuses_sessions_without_a_canonical_request() -> None:
    with pytest.raises(SessionVersionError, match="27 through 30"):
        upgrade_session_document(_fixture("cli-linear-protein.v30.gbdraw-session.json.gz"))


def test_upgrade_keeps_a_settings_only_session_settings_only() -> None:
    source = _fixture("settings-only.v42.json.gz")
    upgrade = upgrade_session_document(source)
    upgraded = upgrade.document

    assert upgrade.warnings == ()
    assert upgraded.version == CURRENT_SESSION_VERSION
    assert SessionDrawing("circular", "Circular", "circular", False) in upgraded.drawings
    assert not any(drawing.has_canonical_request for drawing in upgraded.drawings)
    # The flat draft is split into the mode settings; Circular keeps its values.
    assert "config" not in upgraded.to_dict()
    form = drawing_view(upgraded.to_dict(), "circular")["modes"]["circular"]["config"]["form"]
    assert form and all(source["config"]["form"][key] == value for key, value in form.items())


def test_upgrade_names_and_logs_each_result_it_drops(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # Session 33 saved no feature catalog, so its Results cannot be kept.
    source = _fixture("feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz")
    source["results"].append({**source["results"][0], "name": "second.svg"})
    assert [result["name"] for result in source["results"]] == ["out.svg", "second.svg"]

    with caplog.at_level(logging.WARNING, logger="gbdraw.session"):
        upgrade = upgrade_session_document(source, temporary_directory=tmp_path)

    assert drawing_view(upgrade.document.to_dict(), mode="linear")["results"] == []
    assert upgrade.warnings == (
        f"Upgrading Session 33 to {CURRENT_SESSION_VERSION} dropped the linear drawing's "
        "2 Result(s) 'out.svg', 'second.svg': Session 33 saved no feature catalog "
        "for them. Render the drawing and save the Session to write new Results.",
    )
    assert [record.getMessage() for record in caplog.records] == [
        f"WARNING: {warning}" for warning in upgrade.warnings
    ]
    assert all(record.levelno == logging.WARNING for record in caplog.records)


def test_upgrade_without_results_to_drop_warns_nothing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    source = _fixture("feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz")
    source["results"] = []

    with caplog.at_level(logging.WARNING, logger="gbdraw.session"):
        upgrade = upgrade_session_document(source, temporary_directory=tmp_path)

    assert upgrade.warnings == ()
    assert not caplog.records


def _render_svgs(document: SessionDocument, directory: Path, drawing: str | None = None) -> list[bytes]:
    with materialize_session(document, output_directory=directory) as materialized:
        result = render_session(materialized, drawing=drawing)
    return [path.read_bytes() for path in result.output_paths]


_UPGRADE_FIXTURES = (
    pytest.param("feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz", id="v33-linear-crop"),
    pytest.param("q-frame-main-linear-reverse.v42.gbdraw-session.json.gz", id="v42-reversed-frame"),
    pytest.param("feature-placements-linear.v44.gbdraw-session.json.gz", id="v44-placements"),
)
_SLOW_UPGRADE_FIXTURES = tuple(
    pytest.param(name, id=name.split(".gbdraw-session")[0], marks=pytest.mark.slow)
    for name in (
        "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz",
        "BGC0000708-BGC0000713.v39.gbdraw-session.json.gz",
        "BGC0000708-BGC0000713.v40-schema5.json",
        "HmmtDNA_ATskew.v40-schema5.json",
        "HmmtDNA_basic_circular.issue-469.json.gz",
        "HmmtDNA_basic_circular.v44-schema7.json.gz",
        "HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz",
        "composite-circular-three-files.v44-schema8.gbdraw-session.json.gz",
        "conservation-fasta.v39.gbdraw-session.json.gz",
        "feature-edits-circular-copies.v44.gbdraw-session.json.gz",
        "feature-edits-crop-rc.v44.gbdraw-session.json.gz",
        "feature-placements-circular.v44.gbdraw-session.json.gz",
        "lambda_basic_linear.v40-schema5.json",
        "lambda_basic_linear.v44-schema8.gbdraw-session.json.gz",
        "q-frame-main-web-losatn.v42.gbdraw-session.json.gz",
        "q-frame-main-web-upload.v42.gbdraw-session.json.gz",
        "rendered-v27.v40-schema6.json.gz",
        "se06-main-linear-blast-cli.v42.gbdraw-session.json.gz",
        "se08-main-linear-cli.v42.gbdraw-session.json.gz",
        "selected-feature-annotations.v44.gbdraw-session.json.gz",
        "single.v41-bindings1.json",
        "synthetic_conservation.gbdraw-session.json.gz",
        "test_linear_cli_sidecar_reuses0.v40-schema6.json.gz",
    )
)


@pytest.mark.parametrize("fixture", (*_UPGRADE_FIXTURES, *_SLOW_UPGRADE_FIXTURES))
def test_an_upgraded_session_renders_its_drawing_as_before(fixture: str, tmp_path: Path) -> None:
    source = load_session_document(SESSION_FIXTURES / fixture)
    upgrade = upgrade_session_document(source, temporary_directory=tmp_path / "materialized")
    upgraded = upgrade.document

    assert source.version < upgraded.version == CURRENT_SESSION_VERSION
    # The move also lists the other mode as an undrawn drawing when its settings
    # differ from the defaults (PD-OI-090), as the Web's Load does.
    assert set(source.drawings) <= set(upgraded.drawings)
    assert not any(drawing.has_canonical_request for drawing in set(upgraded.drawings) - set(source.drawings))
    catalog = source.to_dict().get("editorState", {}).get("featureCatalog")
    # Sessions 31-39 saved no catalog: their Results wait for the next render,
    # and the upgrade names each one it drops.
    source_results = source.to_dict()["results"]
    expected_results = source_results if catalog else []
    [source_drawing] = source.drawings
    assert drawing_view(upgraded.to_dict(), source_drawing.id)["results"] == expected_results
    if catalog or not source_results:
        assert upgrade.warnings == ()
    else:
        (warning,) = upgrade.warnings
        assert all(repr(result["name"]) in warning for result in source_results)
    assert _render_svgs(upgraded, tmp_path / "upgraded", source_drawing.id) == _render_svgs(source, tmp_path / "source")


def test_an_upgraded_v33_label_override_with_empty_text_renders_as_before(tmp_path: Path) -> None:
    source = load_session_document(SESSION_FIXTURES / "feature-edits-circular.v33.gbdraw-session.json.gz")
    upgraded = upgrade_session_document(source).document

    assert _render_svgs(upgraded, tmp_path / "upgraded", "circular") == _render_svgs(source, tmp_path / "source")


def test_a_view_keeps_the_key_order_of_its_source() -> None:
    """OV-598: the keys follow the Session they were read from, not a set's iteration order.

    A set of strings iterates in an order that depends on ``PYTHONHASHSEED``;
    one order and its reverse cannot both match it.
    """

    fields = ["renderRequest", "results", "runMetadata", "cliInvocation"]
    editor = ["featureCatalog", "alignmentResetReceipt", "originalSvgStroke"]
    ui = ["selectedResultIndex", "generatedLegendPosition", "appliedPaletteName", "appliedPaletteColors"]
    for reverse in (False, True):
        order = [names[::-1] if reverse else names for names in (fields, editor, ui)]
        source = {
            **{key: {} for key in order[0]},
            "editorState": {**{key: None for key in order[1]}, "legend": {}},
            "ui": {key: 0 for key in order[2]},
        }
        view = _view({}, source)
        assert [key for key in view if key in fields] == order[0]
        assert [key for key in view["editorState"] if key in editor] == order[1]
        assert list(view["ui"]) == order[2]
