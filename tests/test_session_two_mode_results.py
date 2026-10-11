"""Session 46 keeps a Result set of each diagram mode (E1, ``otherModeResult``).

The Web app keeps one Result per mode and saves both: the shown mode's set at
the top level and the other in ``otherModeResult``. The validator admits the
field only beside a committed request of the other mode; the public API names
each set as a drawing (``circular``, ``linear``); the CLI subcommand renders the
drawing of its own mode and a re-save replaces only that drawing's set.
"""

from __future__ import annotations

import copy
import gzip
import json
import re
from pathlib import Path
from typing import Literal, cast

import pytest

from gbdraw.api import (
    CircularDiagramRequest,
    LinearDiagramRequest,
    SessionDrawing,
    SessionDrawingSelectionError,
    load_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.session import session_drawing_artifacts
from gbdraw.cli_utils import session as cli_session
from gbdraw.exceptions import ValidationError
from gbdraw.session_io import load_session, validate_session, write_session_json
from gbdraw.session_resources import canonical_resource_ids
from tests.utils.two_mode_session import two_mode_session


def _feature_ids(svg: str) -> set[str]:
    return set(re.findall(r'data-gbdraw-feature-id="([^"]+)"', svg))


def _first_source_bytes(data: dict) -> str:
    return data["resources"][data["renderRequest"]["records"][0]["source"]["resourceId"]]["data"]


def test_session_46_admits_the_other_mode_result_set() -> None:
    session = two_mode_session()
    validate_session(session)
    assert set(session["otherModeResult"]) <= {
        "renderRequest", "results", "editorState", "ui", "runMetadata", "cliInvocation"
    }


@pytest.mark.parametrize(
    ("change", "message"),
    (
        (lambda s: s.update(version=44), "cannot contain otherModeResult"),
        (lambda s: s["otherModeResult"].update(renderRequest=copy.deepcopy(s["renderRequest"])),
         "committed request of the other mode"),
        (lambda s: s["otherModeResult"].update(results=[]), "Result"),
        (lambda s: s["otherModeResult"].update(config={}), "only a committed Result set"),
        (lambda s: s["otherModeResult"]["ui"].update(mode="linear"), "hold only that Result set"),
        (lambda s: s["otherModeResult"]["renderRequest"]["records"][0]["source"].update(
            resourceId="record-9-genbank"), "missing resource\\(s\\): record-9-genbank"),
    ),
)
def test_session_46_rejects_an_unusable_other_mode_result_set(change, message) -> None:
    session = two_mode_session()
    change(session)
    with pytest.raises(ValidationError, match=message):
        validate_session(session)


def test_settings_only_session_cannot_hold_another_mode_result_set() -> None:
    session = two_mode_session()
    session.update(renderRequest=None, results=[])
    session["editorState"]["featureCatalog"] = None
    with pytest.raises(ValidationError, match="committed request of the other mode"):
        validate_session(session)


def test_session_document_names_each_mode_as_a_drawing(tmp_path: Path) -> None:
    document = load_session_document(two_mode_session())
    assert document.drawings == (
        SessionDrawing("circular", "Circular", "circular", True),
        SessionDrawing("linear", "Linear", "linear", True),
    )
    assert document.active_drawing_id == "circular"
    with pytest.raises(SessionDrawingSelectionError, match="circular .*linear"):
        document.mode
    with pytest.raises(SessionDrawingSelectionError):
        document.has_canonical_request
    assert document.drawing("linear") == document.drawing("Linear") == document.drawing(mode="linear")
    with pytest.raises(SessionDrawingSelectionError, match="no drawing 'linear-2'.*circular .*linear"):
        document.drawing("linear-2")
    with pytest.raises(SessionDrawingSelectionError, match="is a circular drawing, not linear"):
        document.drawing("circular", mode="linear")
    linear = session_drawing_artifacts(document, "linear").fields
    assert linear["renderRequest"]["mode"] == "linear"
    assert "otherModeResult" not in linear
    # Each drawing reads its mode's slice; the slices stay where they are.
    assert linear["modes"] == document.to_dict()["modes"]
    assert linear.get("cliOptions") == document.to_dict().get("cliOptions")
    assert session_drawing_artifacts(document, "circular").fields["renderRequest"] == document.to_dict()["renderRequest"]
    with materialize_session(document, output_directory=tmp_path) as materialized:
        with pytest.raises(SessionDrawingSelectionError):
            session_to_request(materialized)
        assert isinstance(session_to_request(materialized, drawing="linear"), LinearDiagramRequest)
        assert isinstance(session_to_request(materialized, drawing="circular"), CircularDiagramRequest)


@pytest.mark.parametrize("mode", ("circular", "linear"))
def test_cli_session_replay_renders_the_set_of_its_own_mode(tmp_path: Path, mode: str) -> None:
    source = tmp_path / "two.gbdraw-session.json"
    write_session_json(source, two_mode_session())
    output = tmp_path / f"replayed-{mode}"
    assert cli_session.render_canonical_session_if_present(
        load_session(str(source)),
        mode=cast(Literal["circular", "linear"], mode),
        output_override=str(output),
        format_override="svg",
        save_session=False,
        session_output=None,
    )
    session = two_mode_session()
    saved = (session if mode == "circular" else session["otherModeResult"])["results"][0]["content"]
    rendered = _feature_ids(output.with_suffix(".svg").read_text(encoding="utf-8"))
    assert rendered and rendered == _feature_ids(saved)


def test_cli_session_resave_keeps_the_other_mode_result_set(tmp_path: Path) -> None:
    source = tmp_path / "two.gbdraw-session.json"
    session = two_mode_session()
    write_session_json(source, session)
    resaved = tmp_path / "resaved.gbdraw-session.json"
    assert cli_session.render_canonical_session_if_present(
        load_session(str(source)),
        mode="linear",
        output_override=str(tmp_path / "replayed"),
        format_override="svg",
        save_session=False,
        session_output=str(resaved),
    )
    document = load_session_document(resaved)
    # The re-save replaces the Linear drawing's set where it lives and keeps
    # the drawing order and the shown mode (``ui.mode``).
    assert [drawing.id for drawing in document.drawings] == ["circular", "linear"]
    assert document.to_dict()["ui"]["mode"] == session["ui"]["mode"] == "circular"
    circular = session_drawing_artifacts(document, "circular").fields
    assert circular["renderRequest"] == session["renderRequest"]
    assert circular["results"] == session["results"]
    assert circular["editorState"]["featureCatalog"] == session["editorState"]["featureCatalog"]
    # The kept Circular set is unchanged, and every resource it names keeps
    # its ID and bytes in the rebuilt table.
    resaved_resources = document.to_dict()["resources"]
    for resource_id in canonical_resource_ids(session["renderRequest"]):
        assert resaved_resources[resource_id]["data"] == session["resources"][resource_id]["data"]
    # The re-rendered Linear set names its own bytes.
    linear = session_drawing_artifacts(document, "linear").fields
    assert linear["renderRequest"]["mode"] == "linear"
    assert _first_source_bytes(dict(linear)) == _first_source_bytes({
        **session, "renderRequest": session["otherModeResult"]["renderRequest"]
    })
    with materialize_session(document, output_directory=tmp_path) as materialized:
        assert isinstance(session_to_request(materialized, drawing="circular"), CircularDiagramRequest)


SESSION_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "sessions"


@pytest.mark.parametrize(
    "fixture",
    (
        "feature-edits-circular.v33.gbdraw-session.json.gz",
        "BGC0000708-BGC0000713.v39.gbdraw-session.json.gz",
        "conservation-fasta.v39.gbdraw-session.json.gz",
    ),
)
def test_sessions_before_46_reject_another_mode_result_set(fixture: str) -> None:
    session = json.loads(gzip.decompress((SESSION_FIXTURES / fixture).read_bytes()))
    session["otherModeResult"] = two_mode_session()["otherModeResult"]
    with pytest.raises(ValidationError, match="cannot contain otherModeResult"):
        validate_session(session)


@pytest.mark.parametrize(
    ("change", "message"),
    (
        (lambda s: s["otherModeResult"].update(ui=None), "hold only that Result set"),
        (lambda s: s["otherModeResult"]["editorState"].update(legend=None), "hold only that Result set"),
        (lambda s: s["otherModeResult"].update(runMetadata=None), "runMetadata holds only"),
        (lambda s: s["otherModeResult"].update(runMetadata={"annotationWarnings": [{"resultName": "nope", "bogus": 1}]}),
         "Annotation warnings"),
        (lambda s: s["otherModeResult"].update(runMetadata={"comparisonWarnings": {}}), "must be an array"),
    ),
)
def test_the_other_set_has_the_web_readers_rules(change, message) -> None:
    session = two_mode_session()
    change(session)
    with pytest.raises(ValidationError, match=message):
        validate_session(session)


def test_the_other_set_field_lists_match_the_web_reader() -> None:
    """One rule in both validators (services/session-authority.js)."""
    source = (Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "js" / "services" / "session-authority.js").read_text(
        encoding="utf-8"
    )

    def js_set(name: str) -> set[str]:
        body = re.search(rf"const {name} = new Set\(\[(.*?)\]\);", source, re.S)
        assert body, name
        return set(re.findall(r"'([^']+)'", body.group(1)))

    import gbdraw.session_io as session_io

    assert js_set("OTHER_MODE_RESULT_FIELDS") == set(session_io.OTHER_MODE_RESULT_FIELDS)
    assert js_set("OTHER_MODE_RESULT_UI_FIELDS") == set(session_io.OTHER_MODE_RESULT_UI_FIELDS)
    assert js_set("OTHER_MODE_RESULT_EDITOR_FIELDS") == set(session_io.OTHER_MODE_RESULT_EDITOR_FIELDS)
    assert js_set("OTHER_MODE_RESULT_LEGEND_FIELDS") == set(session_io.OTHER_MODE_RESULT_LEGEND_FIELDS)
    assert js_set("OTHER_MODE_RESULT_RUN_METADATA_FIELDS") == set(session_io.OTHER_MODE_RESULT_RUN_METADATA_FIELDS)


def test_a_drawing_view_lacks_the_per_set_fields_its_set_lacks() -> None:
    session = two_mode_session()
    session["ui"].update(appliedPaletteName="tableau", appliedPaletteColors={"CDS": "#111111"}, selectedResultIndex=0)
    session["editorState"]["alignmentResetReceipt"] = None
    session["editorState"].setdefault("legend", {})["originalOrder"] = ["CDS", "repeat_region"]
    other = session["otherModeResult"]
    other["ui"] = {}
    other["editorState"].pop("legend", None)
    other["editorState"].pop("originalSvgStroke", None)
    document = load_session_document(session)
    linear = session_drawing_artifacts(document, "linear").fields
    assert "appliedPaletteName" not in linear["ui"]
    assert "appliedPaletteColors" not in linear["ui"]
    assert "originalOrder" not in linear["editorState"]["legend"]
    assert "originalSvgStroke" not in linear["editorState"]
    circular = session_drawing_artifacts(document, "circular").fields
    assert circular["ui"]["appliedPaletteName"] == "tableau"
    assert circular["editorState"]["legend"]["originalOrder"] == ["CDS", "repeat_region"]
    # The shared draft and preferences stay.
    assert linear["ui"]["mode"] == session["ui"]["mode"]
    assert linear["modes"] == session["modes"]
    assert linear.get("cliOptions") == session.get("cliOptions")
