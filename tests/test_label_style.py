"""The feature label style settings (LP-D98, LP-D100): per-mode defaults, the leader width rule, refusal of
bad values on every surface, and the retired ``labels.stroke_width``."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import Any, Literal

import pytest

import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli
from gbdraw.api import (
    CircularDiagramOptions,
    LinearDiagramOptions,
    load_session_document,
    materialize_session,
    render_session,
)
from gbdraw.api.config import apply_config_overrides, load_default_config
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig, LinearRenderProfile
from gbdraw.config.modify import validate_config_overrides
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import GbdrawError, ValidationError
from gbdraw.labels.label_style import (
    BOX_OPACITY,
    DEFAULT_LABEL_LINE,
    DEFAULT_LABEL_PROTECTION,
    FeatureLabelStyle,
    leader_width_px,
)
from tests.utils.render_stage_equivalence import EquivalenceHarness

REPO_ROOT = Path(__file__).resolve().parents[1]
LAMBDA = REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"
SESSIONS = REPO_ROOT / "tests" / "fixtures" / "sessions"
WEB_OVERRIDES = SESSIONS / "lambda_basic_linear.v44-schema8.gbdraw-session.json.gz"
CLI_CONFIG = SESSIONS / "cli-linear-tables.v44.gbdraw-session.json.gz"


def _session(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def _python_svgs(session: dict[str, Any], directory: Path) -> list[str]:
    with materialize_session(load_session_document(session), output_directory=directory) as materialized:
        render_session(materialized)
    return [path.read_text(encoding="utf-8") for path in sorted(directory.glob("*.svg"))]


def _leader_widths(svg: str) -> set[str]:
    return {
        width
        for tag in re.findall(r"<line [^>]*data-label-feature-id=[^>]*>", svg)
        for width in re.findall(r'stroke-width="([^"]*)"', tag)
    }


def test_each_mode_has_its_packaged_label_style() -> None:
    # Linear: an underline on a white box at 0.8; Circular: neither (a halo can be chosen, LP-D110).
    config = GbdrawConfig.from_dict(load_default_config())
    assert LinearRenderProfile(config).label_style == FeatureLabelStyle("under", "box")
    assert CircularRenderProfile(config).label_style == FeatureLabelStyle("none", "none")
    assert BOX_OPACITY == 0.8
    packaged = load_config_toml("gbdraw.data", "config.toml")["labels"]
    for mode in ("linear", "circular"):
        assert (packaged[mode]["line"], packaged[mode]["protection"]) == (
            DEFAULT_LABEL_LINE[mode],
            DEFAULT_LABEL_PROTECTION[mode],
        )
        # A stored config written before the settings existed takes the same defaults.
        stored = load_default_config()
        del stored["labels"][mode]["line"], stored["labels"][mode]["protection"]
        labels = getattr(GbdrawConfig.from_dict(stored).labels, mode)
        assert (labels.line, labels.protection) == (DEFAULT_LABEL_LINE[mode], DEFAULT_LABEL_PROTECTION[mode])


@pytest.mark.parametrize(("font_size", "width"), [(5.0, 1.0), (8.0, 1.0), (14.0, 1.0), (24.0, 1.68)])
def test_leaders_and_underlines_are_007_em_wide_and_at_least_1_px(font_size: float, width: float) -> None:
    assert leader_width_px(font_size) == pytest.approx(width)


def _linear_label_with_leader(font_size: float, x: float, leader_end_y: float) -> dict[str, object]:
    """A prepared Linear label above its feature, its leader running straight down from y 0 to ``leader_end_y``."""

    return {
        "label_text": "orf",
        "font_family": "'Liberation Sans', sans-serif",
        "font_size": font_size,
        "middle_x": x,
        "middle_y": -2.0 * font_size,
        "width_px": 2.0 * font_size,
        "height_px": font_size,
        "rotation_deg": 0.0,
        "text_anchor": "middle",
        "is_embedded": False,
        "leader_line": True,
        "leader_start_x": x,
        "leader_start_y": 0.0,
        "leader_end_x": x,
        "leader_end_y": leader_end_y,
    }


def test_each_leader_keeps_its_own_width_clear_in_a_figure_with_two_label_sizes() -> None:
    """A 24 px label's leader is 1.68 px wide and a 5 px label's 1 px: the label band, the region-label
    obstacles of both modes and the drawing read each label's own width, not the widest one (review LS-1, 6)."""

    from types import SimpleNamespace

    from gbdraw.diagrams.linear.annotation_labels import linear_feature_label_obstacles
    from gbdraw.layout.linear import measure_linear_label_band

    big, small = _linear_label_with_leader(24.0, 0.0, 30.0), _linear_label_with_leader(5.0, 100.0, 40.0)
    band = measure_linear_label_band([big, small])
    assert band is not None and band.bottom_y == pytest.approx(40.0 + 0.5)

    style = FeatureLabelStyle(line="none", protection="none")
    leaders = {
        obstacle.owner[2]: obstacle.box.max_x - obstacle.box.min_x
        for obstacle in linear_feature_label_obstacles([big, small], SimpleNamespace(slots=()), 0, label_style=style)
        if obstacle.leader
    }
    assert leaders == {0: pytest.approx(1.68), 1: pytest.approx(1.0)}

    # A label 3.6 px from a 1.68 px leader and from a 1 px leader, keeping 3 px: it meets only the wide one.
    from gbdraw.render.drawers.circular.annotations import CircularLabelObstacles, LeaderSegment, _ObstacleIndex

    square = ((0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0))
    obstacles = CircularLabelObstacles(
        segments=(
            LeaderSegment((13.6, -5.0), (13.6, 15.0), 0.5 * leader_width_px(24.0)),
            LeaderSegment((-3.6, -5.0), (-3.6, 15.0), 0.5 * leader_width_px(5.0)),
        ),
        clearance=3.0,
    )
    assert _ObstacleIndex(obstacles).conflicts(SimpleNamespace(parts=(square,))) == 1


def test_only_a_horizontal_label_gets_the_underline_and_its_obstacles_agree() -> None:
    """LP-D112: a feature label drawn turned gets no underline; its leader ends at its letters and what the region
    labels keep off is what ``line = none`` draws. A horizontal label with the same style keeps its underline."""

    from types import SimpleNamespace

    from gbdraw.diagrams.linear.annotation_labels import linear_feature_label_obstacles
    from gbdraw.labels.label_ink import linear_label_ink
    from gbdraw.labels.linear import linear_leader_polyline
    from gbdraw.render.drawers.label_backing import feature_label_drawn_parts

    under, plain = FeatureLabelStyle(line="under", protection="box"), FeatureLabelStyle(line="none", protection="box")
    assert under.draws_underline(turned=False) and not under.draws_underline(turned=True)
    assert not plain.draws_underline(turned=False)
    flat = _linear_label_with_leader(10.0, 0.0, -14.0)
    turned = {**flat, "rotation_deg": -45.0, "text_anchor": "start"}
    flat_ink = linear_label_ink(flat, linear_leader_polyline(flat), under)
    assert flat_ink.underline is not None and flat_ink.underline_box is not None
    # A full turn draws the label horizontal again, so it keeps the underline.
    assert linear_label_ink({**flat, "rotation_deg": -360.0}, linear_leader_polyline(flat), under).underline is not None
    ink = linear_label_ink(turned, leader := linear_leader_polyline(turned), under)
    assert ink.underline is None and ink.underline_box is None
    assert ink == linear_label_ink(turned, leader, plain)
    assert feature_label_drawn_parts(ink, under) == feature_label_drawn_parts(ink, plain)
    plan = SimpleNamespace(slots=())
    assert linear_feature_label_obstacles([turned], plan, 0, label_style=under) == linear_feature_label_obstacles(
        [turned], plan, 0, label_style=plain
    )
    assert linear_feature_label_obstacles([flat], plan, 0, label_style=under) != linear_feature_label_obstacles(
        [flat], plan, 0, label_style=plain
    )


@pytest.mark.parametrize("angle", [90.0, 270.0])
def test_a_circular_radial_label_gets_no_underline_even_where_it_reads_horizontally(angle: float) -> None:
    """LP-D112: a Circular radial label is turned whatever its angle; at 3 and 9 o'clock its rotation is 0 (it
    reads horizontally) and it still gets no underline, while a horizontal (arc) label there gets one."""

    from gbdraw.labels.circular_radial import radial_text_orientation
    from gbdraw.labels.label_ink import circular_label_ink

    under = FeatureLabelStyle(line="under", protection="none")
    rotation, anchor, baseline = radial_text_orientation(angle, is_inner=False)
    point = {"start_x": 100.0, "start_y": 0.0, "middle_x": 90.0, "middle_y": 0.0, "feature_middle_x": 80.0}
    label = {
        "label_text": "orf",
        "font_family": "'Liberation Sans', sans-serif",
        "font_size": 10.0,
        "middle": 1000,
        "is_inner": False,
        **point,
    }
    radial = {**label, "placement": "radial", "rotation_deg": rotation, "text_anchor": anchor, "dominant_baseline": baseline}
    assert rotation == 0.0
    assert circular_label_ink(radial, 4000, under).underline is None
    assert circular_label_ink(label, 4000, under).underline is not None


@pytest.mark.parametrize(
    ("options", "path", "value"),
    [
        (LinearDiagramOptions, "labels.linear.line", "over"),
        (LinearDiagramOptions, "labels.linear.protection", "shadow"),
        (CircularDiagramOptions, "labels.circular.line", "over"),
        (CircularDiagramOptions, "labels.circular.protection", "shadow"),
    ],
)
def test_a_bad_style_value_is_refused_by_python_the_cli_and_a_session(
    options: type, path: str, value: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(ValidationError, match=re.escape(path)) as caught:
        options(config_overrides={path: value})
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "configPath": path}

    cli = linear_cli if path.startswith("labels.linear") else circular_cli
    flag = "--label_line" if path.endswith(".line") else "--label_protection"
    with pytest.raises(SystemExit):
        cli._get_args(["--gbk", str(LAMBDA), flag, value])
    assert "invalid choice" in capsys.readouterr().err

    if path.startswith("labels.linear"):
        session = _session(WEB_OVERRIDES)
        session["renderRequest"]["diagramOptions"]["configOverrides"][path] = value
        with pytest.raises(GbdrawError, match=re.escape(path)):
            _python_svgs(session, tmp_path)


def test_the_retired_stroke_width_override_is_an_unknown_path() -> None:
    for path in ("labels.stroke_width.long", "labels.stroke_width.short"):
        with pytest.raises(ValidationError) as caught:
            validate_config_overrides({path: 0.5})
        assert caught.value.diagnostic["reason"] == "UNKNOWN_CONFIG_PATH"
        with pytest.raises(ValidationError, match="Unknown config override path"):
            LinearDiagramOptions(config_overrides={path: 0.5})


def test_a_stored_config_with_the_old_stroke_width_draws_the_new_label_style(tmp_path: Path) -> None:
    # A CLI Session stores its full config, with the old default width (no Gallery Session stores one since #1217).
    session = _session(CLI_CONFIG)
    stored = session["renderRequest"]["diagramOptions"]["config"]["labels"]
    assert stored["stroke_width"] == {"short": 1.5, "long": 0.5}
    svgs = _python_svgs(session, tmp_path)
    assert svgs
    for svg in svgs:
        widths = _leader_widths(svg)
        assert widths and "0.5" not in widths and "1.5" not in widths, widths
        # Linear draws its default box, not a halo (LP-D106).
        assert "<use" not in svg and f'fill-opacity="{BOX_OPACITY}"' in svg


def test_the_four_settings_reach_the_drawing_from_cli_and_python(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    linear_cli.linear_main(
        ["--gbk", str(LAMBDA), "--show_labels", "--label_line", "none", "--label_protection", "box", "-o", "l", "-f", "svg"]
    )
    circular_cli.circular_main(
        ["--gbk", str(LAMBDA), "--labels", "--label_line", "under", "--label_protection", "none", "-o", "c", "-f", "svg"]
    )
    linear, circular = (tmp_path / "l.svg").read_text(), (tmp_path / "c.svg").read_text()
    assert "<use" not in linear and f'fill-opacity="{BOX_OPACITY}"' in linear
    assert "<use" not in circular and f'fill-opacity="{BOX_OPACITY}"' not in circular
    cfg = apply_config_overrides(
        GbdrawConfig.from_dict(load_default_config()),
        {"labels.linear.line": "none", "labels.circular.protection": "box"},
    )
    assert LinearRenderProfile(cfg).label_style == FeatureLabelStyle("none", "box")
    assert CircularRenderProfile(cfg).label_style == FeatureLabelStyle("none", "box")


@pytest.mark.parametrize(
    ("mode", "argv"),
    [
        ("linear", ["--show_labels", "--separate_strands"]),
        ("circular", ["--labels", "both"]),
    ],
)
def test_a_style_change_rebuilds_its_stages_without_an_undeclared_read(
    mode: Literal["linear", "circular"], argv: list[str], tmp_path: Path
) -> None:
    harness = EquivalenceHarness(mode, tmp_path)
    base = [*argv, "--gbk", str(LAMBDA)]
    harness.warm(base)
    other = {"linear": "halo", "circular": "box"}[mode]
    for extra in (["--label_line", "under" if mode == "circular" else "none"], ["--label_protection", other]):
        warm, _cold = harness.step(base + extra, " ".join(extra))
        assert warm.error is None, warm.error
        assert warm.stages.undeclared_read_total == 0 and warm.stages.uncacheable_total == 0, warm.stages
        # The style is in the keys of the stages that draw it: they are rebuilt, not reused.
        rebuilt = {"record_content"} if mode == "linear" else {"label_leaders_content", "label_text_content"}
        assert rebuilt & set(warm.stages.builds), warm.stages.builds


def test_a_label_gets_its_box_when_it_has_a_leader_or_sits_beside_its_feature() -> None:
    """LP-D121 (OV-1361, LP-D119): ``protection = box`` gives a feature label its white box when the label has
    a leader or is drawn on its feature without one (Linear ``above_feature``), turned or not; a label drawn
    inside its feature has neither and gets no box, and ``protection = none`` draws none."""

    from gbdraw.labels.label_ink import linear_label_ink
    from gbdraw.labels.linear import linear_leader_polyline

    box, none = FeatureLabelStyle(line="under", protection="box"), FeatureLabelStyle(line="under", protection="none")
    for has_leader in (False, True):
        for beside_feature in (False, True):
            assert box.draws_box(has_leader=has_leader, beside_feature=beside_feature) is (has_leader or beside_feature)
            assert not none.draws_box(has_leader=has_leader, beside_feature=beside_feature)

    with_leader = _linear_label_with_leader(10.0, 0.0, -14.0)
    inside = {**with_leader, "is_embedded": True, "leader_line": False, "middle_y": 0.0}
    beside = {**inside, "middle_y": -8.0, "above_feature_place_above": True}
    for rotation in (0.0, -45.0):
        drawn = {
            name: linear_label_ink(turned, linear_leader_polyline(turned), box)
            for name, label in (("leader", with_leader), ("beside", beside), ("inside", inside))
            for turned in [{**label, "rotation_deg": rotation}]
        }
        assert drawn["leader"].draws_box(box)
        assert drawn["beside"].draws_box(box) and drawn["beside"].feature_edge_y == 0.0
        assert not drawn["inside"].draws_box(box) and drawn["inside"].feature_edge_y is None
