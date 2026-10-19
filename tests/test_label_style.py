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
    # Linear: an underline on a white box at 0.8; Circular: a halo (LP-D106).
    config = GbdrawConfig.from_dict(load_default_config())
    assert LinearRenderProfile(config).label_style == FeatureLabelStyle("under", "box")
    assert CircularRenderProfile(config).label_style == FeatureLabelStyle("none", "halo")
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
