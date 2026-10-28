"""Hide the Circular centre definition (PD-OI-109).

``objects.definition.circular.show`` (default true; CLI ``--hide_definition``)
turns every record's centre definition off. A hidden definition draws no
record-definition group, reserves no band for the inside tracks, and raises no
definition fit warning; the plot title is drawn as before.
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble
import gbdraw.labels.circular_fit as circular_fit_module
from gbdraw.api import CircularDiagramOptions, CircularMultiRecordOptions, CircularOutputOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram
from gbdraw.config.modify import canonical_config_override_paths

TEST_INPUTS = Path(__file__).parent / "test_inputs"
HMMT = TEST_INPUTS / "HmmtDNA.gbk"
MERS_COV = TEST_INPUTS / "MERS-CoV.gbk"
LAMBDA = TEST_INPUTS / "NC_001416.gb"
SVG_NS = {"svg": "http://www.w3.org/2000/svg"}
SHOW_PATH = "objects.definition.circular.show"
HIDDEN = {SHOW_PATH: False}


def _record(path: Path = HMMT):
    return next(SeqIO.parse(str(path), "genbank"))


def _role_groups(svg_text: str, role: str) -> list[ET.Element]:
    return ET.fromstring(svg_text).findall(f".//svg:g[@data-gbdraw-role='{role}']", SVG_NS)


def _track_paths(svg_text: str) -> dict[str, str]:
    """The drawn feature, GC content and GC skew paths by group id: their radii and widths."""

    paths: dict[str, str] = {}
    for group in ET.fromstring(svg_text).iter("{http://www.w3.org/2000/svg}g"):
        group_id = str(group.get("id") or "")
        if group_id.startswith(("gc_content", "skew", "gc_skew", "tick")):
            paths[group_id] = "|".join(str(path.get("d")) for path in group.iter("{http://www.w3.org/2000/svg}path"))
    assert paths
    return paths


def test_the_packaged_config_shows_the_definition() -> None:
    assert apply_config_overrides(None, None).objects.definition.circular.show is True
    assert apply_config_overrides(None, HIDDEN).objects.definition.circular.show is False
    assert SHOW_PATH in canonical_config_override_paths()


@pytest.mark.circular
def test_a_hidden_definition_draws_no_definition_for_one_record() -> None:
    shown = build_circular_diagram(_record()).tostring()
    hidden = build_circular_diagram(_record(), options=CircularDiagramOptions(config_overrides=HIDDEN)).tostring()

    assert _role_groups(shown, "record-definition")
    assert _role_groups(hidden, "record-definition") == []


@pytest.mark.circular
def test_a_hidden_definition_draws_no_definition_on_a_multi_record_canvas() -> None:
    records = [_record(), _record(MERS_COV)]
    layout = CircularMultiRecordOptions(multi_record_positions=["#1@1", "#2@1"])
    shown = build_circular_multi_diagram(records, layout=layout).tostring()
    hidden = build_circular_multi_diagram(
        records, options=CircularDiagramOptions(config_overrides=HIDDEN), layout=layout
    ).tostring()

    assert len(_role_groups(shown, "record-definition")) == 2
    assert _role_groups(hidden, "record-definition") == []


@pytest.mark.circular
@pytest.mark.parametrize("multi_record", [False, True], ids=["one record", "multi-record canvas"])
def test_a_hidden_definition_keeps_the_plot_title(multi_record: bool) -> None:
    options = CircularDiagramOptions(
        config_overrides=HIDDEN,
        output=CircularOutputOptions(plot_title_position="top"),
        plot_title="Mitochondrial genome",
    )
    drawing = (
        build_circular_multi_diagram([_record(), _record(MERS_COV)], options=options)
        if multi_record
        else build_circular_diagram(_record(), options=options)
    )
    svg_text = drawing.tostring()
    titles = _role_groups(svg_text, "plot-title")

    assert titles
    assert "Mitochondrial genome" in "".join(titles[0].itertext())
    assert _role_groups(svg_text, "record-definition") == []


@pytest.mark.circular
def test_the_shown_definition_override_draws_the_default_figure() -> None:
    default = build_circular_diagram(_record()).tostring()
    shown = build_circular_diagram(
        _record(), options=CircularDiagramOptions(config_overrides={SHOW_PATH: True})
    ).tostring()

    assert shown == default


@pytest.mark.circular
def test_a_hidden_definition_reserves_no_centre_band(monkeypatch: pytest.MonkeyPatch) -> None:
    """On the Web default tracks, lambda's definition moves or narrows the GC rows; hidden, they are placed as for a
    record whose definition has no band."""

    web_tracks = {
        "canvas.circular.track_type": "tuckin",
        "canvas.strandedness": True,
        "canvas.show_gc": True,
        "canvas.show_skew": True,
    }

    def tracks(overrides: dict[str, object]) -> dict[str, str]:
        options = CircularDiagramOptions(config_overrides={**web_tracks, **overrides})
        return _track_paths(build_circular_diagram(_record(LAMBDA), options=options).tostring())

    shown = tracks({})
    hidden = tracks(HIDDEN)
    monkeypatch.setattr(circular_assemble, "_definition_reserved_radius_px", lambda *args, **kwargs: 0.0)
    no_band = tracks({})

    assert hidden == no_band
    assert hidden != shown


@pytest.mark.circular
def test_cli_hide_definition_draws_no_definition_and_no_fit_warning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """MERS-CoV's definition does not fit at the default font and warns; ``--hide_definition`` draws none and
    warns nothing."""

    from gbdraw import cli

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)

    def draw(name: str, *args: str) -> tuple[str, str]:
        output = tmp_path / name
        monkeypatch.setattr(
            sys, "argv", ["gbdraw", "circular", "--gbk", str(MERS_COV), *args, "-o", str(output), "-f", "svg"]
        )
        cli.main()
        return output.with_suffix(".svg").read_text(), capsys.readouterr().err

    shown_svg, shown_err = draw("shown")
    hidden_svg, hidden_err = draw("hidden", "--hide_definition")

    assert "definition line does not fit" in shown_err
    assert _role_groups(shown_svg, "record-definition")
    assert hidden_err == ""
    assert _role_groups(hidden_svg, "record-definition") == []


def test_cli_parses_hide_definition() -> None:
    from gbdraw.circular import _get_args

    assert _get_args(["--gbk", "dummy.gb"]).hide_definition is False
    assert _get_args(["--gbk", "dummy.gb", "--hide_definition"]).hide_definition is True


@pytest.mark.circular
@pytest.mark.parametrize(
    ("fixture", "definitions"),
    [
        # The Web writer's Session 49 with Show Definition off.
        ("hidden-definition-circular.v49.gbdraw-session.json.gz", 0),
        # A Session 48 has no setting and draws its definition as before.
        ("record-identity-circular-single.v48.gbdraw-session.json.gz", 1),
    ],
)
def test_a_session_replays_its_definition_setting(
    fixture: str, definitions: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gbdraw import cli

    output = tmp_path / "replay"
    session = Path(__file__).parent / "fixtures" / "sessions" / fixture
    monkeypatch.setattr(
        sys, "argv", ["gbdraw", "circular", "--session", str(session), "-o", str(output), "-f", "svg"]
    )
    cli.main()

    assert len(_role_groups(output.with_suffix(".svg").read_text(), "record-definition")) == definitions
