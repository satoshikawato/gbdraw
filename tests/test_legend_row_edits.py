"""Python draws and fits Legend row edits and Canvas padding (issue #1427, PR-D).

The config leaves ``objects.legends.row_names``, ``hidden_rows`` and
``row_order`` rename, hide and order Legend rows by key; ``canvas.padding``
adds room per side outside the edge padding. ``--legend_rows FILE`` and
``--canvas_padding TOP,RIGHT,BOTTOM,LEFT`` set them from the CLI, a saved
Session keeps them, and a requested figure size fits the edited Legend and the
padding. How the drawn rows match the Web's executor is pinned by
``tests/fixtures/legend_row_edit_vectors.json``.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path

import pytest

import gbdraw
from gbdraw.api.config import apply_config_overrides
from gbdraw.config.models.canvas import CANVAS_PADDING_PATH
from gbdraw.config.models.objects import LEGEND_HIDDEN_ROWS_PATH, LEGEND_ROW_NAMES_PATH, LEGEND_ROW_ORDER_PATH
from gbdraw.exceptions import ValidationError
from gbdraw.render.composition import applied_composition_plan

INPUTS = Path(__file__).resolve().parent / "test_inputs"
MITO = INPUTS / "HmmtDNA.gbk"
LONG_NAME = "Protein-coding genes of the human mitochondrial genome"
RENAMED = "Protein-coding genes"
PADDING = (10.0, 20.0, 30.0, 40.0)


def _mito() -> list:
    return gbdraw.read_genbank(MITO)


def _draw(mode: str, **options) -> gbdraw.Diagram:
    if mode == "circular":
        return gbdraw.draw_circular(_mito(), options=gbdraw.CircularOptions(**options))
    return gbdraw.draw_linear(_mito(), options=gbdraw.LinearOptions(**options))


def _rows(svg: str) -> list[tuple[str, str | None, str]]:
    """The feature Legend rows of the first feature group: key, display, and text."""

    rows = []
    for match in re.finditer(r'<g ([^>]*data-legend-key="[^"]*"[^>]*)>(.*?)</g>', svg):
        attributes = dict(re.findall(r'([\w:-]+)="([^"]*)"', match[1]))
        text = re.search(r"<text[^>]*>([^<]*)</text>", match[2])
        rows.append((html.unescape(attributes["data-legend-key"]), attributes.get("display"),
                     html.unescape(text[1]) if text else ""))
    return rows


def _composition_metadata(svg: str) -> dict:
    head = svg[svg.index("<svg") : svg.index(">", svg.index("<svg"))]
    found = re.search(r'data-gbdraw-composition="([^"]+)"', head)
    assert found
    return json.loads(html.unescape(found[1]))


def _legend_bounds(diagram: gbdraw.Diagram):
    plan = applied_composition_plan(diagram._drawing)
    return plan, next(placement.final_bounds for placement in plan.placements if placement.role == "legend")


# --- the leaves ---


def test_the_leaves_are_empty_by_default_and_typed_when_set() -> None:
    config = apply_config_overrides(None, None)
    assert (config.objects.legends.row_names, config.objects.legends.hidden_rows, config.objects.legends.row_order) == (
        (), (), ())
    assert config.canvas.padding == (0.0, 0.0, 0.0, 0.0)
    edited = apply_config_overrides(None, {
        LEGEND_ROW_NAMES_PATH: [["CDS", "Genes"]], LEGEND_HIDDEN_ROWS_PATH: ["tRNA"],
        LEGEND_ROW_ORDER_PATH: ["rRNA", "CDS"], CANVAS_PADDING_PATH: [1, 2.5, 0, 4],
    })
    assert edited.objects.legends.row_names == (("CDS", "Genes"),)
    assert edited.objects.legends.hidden_rows == ("tRNA",)
    assert edited.objects.legends.row_order == ("rRNA", "CDS")
    assert edited.canvas.padding == (1.0, 2.5, 0.0, 4.0)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (LEGEND_ROW_NAMES_PATH, "CDS"),
        (LEGEND_ROW_NAMES_PATH, [["CDS"]]),
        (LEGEND_ROW_NAMES_PATH, [["CDS", " "]]),
        (LEGEND_HIDDEN_ROWS_PATH, "tRNA"),
        (LEGEND_HIDDEN_ROWS_PATH, [""]),
        (LEGEND_ROW_ORDER_PATH, [1]),
        (CANVAS_PADDING_PATH, [1, 2, 3]),
        (CANVAS_PADDING_PATH, [1, 2, 3, -1]),
        (CANVAS_PADDING_PATH, [1, 2, 3, float("inf")]),
        (CANVAS_PADDING_PATH, [True, 2, 3, 4]),
        (CANVAS_PADDING_PATH, "1,2,3,4"),
    ],
)
def test_a_bad_leaf_is_refused_with_its_config_path(path: str, value: object) -> None:
    with pytest.raises(ValidationError) as error:
        apply_config_overrides(None, {path: value})
    assert error.value.diagnostic == {"code": "INPUT_INVALID", "configPath": path}


# --- the drawn rows ---


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_each_edit_kind_is_drawn(mode: str) -> None:
    unedited = _rows(_draw(mode).to_svg())
    keys = [key for key, _display, _text in unedited]
    assert {"CDS", "tRNA"} <= set(keys)
    renamed = _rows(_draw(mode, config_overrides={LEGEND_ROW_NAMES_PATH: [["CDS", "Genes"]]}).to_svg())
    assert ("Genes", None, "Genes") in renamed and all(key != "CDS" for key, _d, _t in renamed)
    hidden = _rows(_draw(mode, config_overrides={LEGEND_HIDDEN_ROWS_PATH: ["tRNA"]}).to_svg())
    assert ("tRNA", "none", "tRNA") in hidden
    order = list(reversed(dict.fromkeys(keys)))
    ordered = _rows(_draw(mode, config_overrides={LEGEND_ROW_ORDER_PATH: order}).to_svg())
    assert list(dict.fromkeys(key for key, _d, _t in ordered)) == order


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_requests_without_the_leaves_draw_as_before(mode: str) -> None:
    """A request with the empty leaves is the request without them: the same bytes."""

    plain = _draw(mode).to_svg()
    empty = _draw(mode, config_overrides={LEGEND_ROW_NAMES_PATH: [], LEGEND_HIDDEN_ROWS_PATH: [],
                                          LEGEND_ROW_ORDER_PATH: [], CANVAS_PADDING_PATH: [0, 0, 0, 0]}).to_svg()
    assert empty == plain
    assert "canvasPadding" not in _composition_metadata(plain)["spacing"]


# --- the fit and the composition ---


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_auto_canvas_padding_grows_the_canvas_by_the_padding(mode: str) -> None:
    auto = applied_composition_plan(_draw(mode)._drawing)
    padded_diagram = _draw(mode, config_overrides={CANVAS_PADDING_PATH: list(PADDING)})
    padded = applied_composition_plan(padded_diagram._drawing)
    top, right, bottom, left = PADDING
    assert padded.canvas_bounds.width == pytest.approx(auto.canvas_bounds.width + left + right)
    assert padded.canvas_bounds.height == pytest.approx(auto.canvas_bounds.height + top + bottom)
    assert _composition_metadata(padded_diagram.to_svg())["spacing"]["canvasPadding"] == {
        "top": top, "right": right, "bottom": bottom, "left": left}


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_sized_figure_fits_a_longer_renamed_row_and_the_padding(mode: str) -> None:
    width = 250 * 96 / 25.4
    plain_plan, plain_legend = _legend_bounds(_draw(mode, figure_width="250mm"))
    diagram = _draw(mode, figure_width="250mm", config_overrides={
        LEGEND_ROW_NAMES_PATH: [["CDS", RENAMED]], CANVAS_PADDING_PATH: list(PADDING)})
    plan, legend = _legend_bounds(diagram)
    _top, right, _bottom, left = PADDING
    assert plan.canvas_bounds.width == pytest.approx(width, abs=0.01)
    assert plain_plan.canvas_bounds.width == pytest.approx(width, abs=0.01)
    # The fit measured the renamed row: the Legend is wider, and it and the padding stay inside the canvas.
    assert legend.width > plain_legend.width
    assert plan.content_bounds.width <= width + 1e-6
    assert plan.canvas_bounds.min_x - 1e-6 <= plan.content_bounds.min_x
    assert legend.min_x >= plan.content_bounds.min_x + left - 1e-6
    assert legend.max_x <= plan.content_bounds.max_x - right + 1e-6
    assert RENAMED in diagram.to_svg()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_refusal_counts_the_renamed_row_and_the_padding(mode: str) -> None:
    def refusal(overrides: dict[str, object]) -> str:
        with pytest.raises(ValidationError) as error:
            _draw(mode, figure_width="60mm", config_overrides=overrides)
        return str(error.value)

    def part(message: str, name: str) -> float:
        found = re.search(rf"{name} ([0-9.]+) mm", message)
        assert found, message
        return float(found[1])

    plain = refusal({})
    edited = refusal({LEGEND_ROW_NAMES_PATH: [["CDS", LONG_NAME]], CANVAS_PADDING_PATH: list(PADDING)})
    assert part(edited, "legend") > part(plain, "legend")
    # The padding part grows by the left and right canvas padding (60 px).
    assert part(edited, "padding") == pytest.approx(part(plain, "padding") + 60 * 25.4 / 96, abs=0.1)


# --- the CLI and the Session ---


def _rows_table(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "rows.tsv"
    path.write_text(text, encoding="utf-8")
    return path


ROWS_TSV = "key\tname\thidden\tposition\nCDS\tGenes\tfalse\t2\ntRNA\t\ttrue\t3\nrRNA\t\t\t1\n"


def _gbdraw(monkeypatch, capsys, *args: str) -> tuple[int, str]:
    """Run the ``gbdraw`` command in this process: its exit status and its standard error."""

    from gbdraw import cli

    monkeypatch.setattr(sys, "argv", ["gbdraw", *args])
    try:
        cli.main()
        status = 0
    except SystemExit as exit_:
        status = int(exit_.code or 0)
    return status, capsys.readouterr().err


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_cli_options_draw_the_edits_and_a_saved_session_replays_them(
    mode: str, tmp_path: Path, monkeypatch, capsys
) -> None:
    rows = _rows_table(tmp_path, ROWS_TSV)
    out = tmp_path / "out"
    session = tmp_path / "out.gbdraw-session.json"
    status, err = _gbdraw(monkeypatch, capsys, mode, "--gbk", str(MITO), "--legend_rows", str(rows),
                          "--canvas_padding", "10,20,30,40", "-o", str(out), "-f", "svg", "--save_session",
                          "--session_output", str(session))
    assert status == 0, err
    drawn = (tmp_path / "out.svg").read_text(encoding="utf-8")
    assert [key for key, _d, _t in _rows(drawn)][:3] == ["rRNA", "Genes", "tRNA"]
    assert ("tRNA", "none", "tRNA") in _rows(drawn)
    assert _composition_metadata(drawn)["spacing"]["canvasPadding"] == {
        "top": 10.0, "right": 20.0, "bottom": 30.0, "left": 40.0}
    replay = tmp_path / "replay"
    status, err = _gbdraw(monkeypatch, capsys, "render", "--session", str(session), "-o", str(replay), "-f", "svg")
    assert status == 0, err
    assert (tmp_path / "replay.svg").read_text(encoding="utf-8") == drawn


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("CDS\tGenes\n", "the first line must be the header"),
        ("key\tname\thidden\tposition\nCDS\t\tmaybe\t\n", "hidden must be true or false"),
        ("key\tname\thidden\tposition\nCDS\t\t\t0\n", "position must be a whole number from 1"),
        ("key\tname\thidden\tposition\nCDS\t\t\t1\ntRNA\t\t\t\n", "give a position on every row or on none"),
        ("key\tname\thidden\tposition\nCDS\t\t\t1\ntRNA\t\t\t1\n", "each position must be on one row"),
        ("key\tname\thidden\tposition\nCDS\tA\t\t\nCDS\tB\t\t\n", "each key must be on one row"),
    ],
)
def test_the_cli_refuses_a_bad_rows_table(text: str, message: str, tmp_path: Path, monkeypatch, capsys) -> None:
    rows = _rows_table(tmp_path, text)
    status, err = _gbdraw(monkeypatch, capsys, "linear", "--gbk", str(MITO), "--legend_rows", str(rows),
                          "-o", str(tmp_path / "x"))
    assert status == 2
    assert "--legend_rows" in err and message in err


@pytest.mark.parametrize("value", ["1,2,3", "1,2,3,-4", "a,b,c,d"])
def test_the_cli_refuses_a_bad_padding(value: str, tmp_path: Path, monkeypatch, capsys) -> None:
    status, err = _gbdraw(monkeypatch, capsys, "circular", "--gbk", str(MITO), "--canvas_padding", value,
                          "-o", str(tmp_path / "x"))
    assert status == 2
    assert "--canvas_padding" in err


def test_a_saved_session_holds_the_new_config_fields_only_when_set(tmp_path: Path, monkeypatch, capsys) -> None:
    """A Session without the edits is the one written before the fields existed (its full config lacks them)."""

    def saved_config(*args: str) -> dict:
        session = tmp_path / f"s{len(args)}.json"
        status, err = _gbdraw(monkeypatch, capsys, "linear", "--gbk", str(MITO), *args, "-o", str(tmp_path / "o"),
                              "-f", "svg", "--save_session", "--session_output", str(session), "--overwrite")
        assert status == 0, err
        document = json.loads(session.read_text(encoding="utf-8"))
        return document["drawings"][0]["renderRequest"]["diagramOptions"]["config"]

    plain = saved_config()
    assert "padding" not in plain["canvas"]
    assert not {"row_names", "hidden_rows", "row_order"} & set(plain["objects"]["legends"])
    edited = saved_config("--legend_rows", str(_rows_table(tmp_path, ROWS_TSV)), "--canvas_padding", "0,0,5,0")
    assert edited["canvas"]["padding"] == [0.0, 0.0, 5.0, 0.0]
    assert edited["objects"]["legends"]["row_names"] == [["CDS", "Genes"]]
    assert edited["objects"]["legends"]["hidden_rows"] == ["tRNA"]
    assert edited["objects"]["legends"]["row_order"] == ["rRNA", "CDS", "tRNA"]
