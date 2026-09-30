"""GE-03: the Circular definition interval follows its font unless set."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

from gbdraw import cli
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import _resolve_diagram_options_config
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.api.request_render import render_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.config.models.objects import circular_definition_interval_for_font
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request

REPO_ROOT = Path(__file__).resolve().parents[1]
RECORD = REPO_ROOT / "tests" / "fixtures" / "regex_rules.gb"
FONT = "objects.definition.circular.font_size"
INTERVAL = "objects.definition.circular.interval"


@pytest.mark.parametrize(("font", "interval"), [(18, 20), (28, 30), (30, 32), (12.6, 14)])
def test_interval_rule_is_owned_once(font, interval):
    assert circular_definition_interval_for_font(font) == interval


@pytest.mark.parametrize(
    ("overrides", "expected_interval"),
    [
        ({FONT: 30}, 32),
        ({FONT: 30, INTERVAL: 20}, 20),
        ({FONT: 18}, 20),
        ({}, 20),
    ],
)
def test_applied_overrides_derive_interval_only_when_unset(overrides, expected_interval):
    config = apply_config_overrides(None, overrides)
    assert config.objects.definition.circular.interval == expected_interval
    options = CircularDiagramOptions(config_overrides=overrides)
    resolved = _resolve_diagram_options_config(options)
    assert resolved.objects.definition.circular.interval == expected_interval


def _decode(tmp_path: Path, overrides: dict):
    encoded = encode_canonical_request(
        CircularDiagramRequest(records=(RecordInput(source=GenBankInputSource(RECORD)),))
    )
    payload = copy.deepcopy(encoded.payload)
    payload["diagramOptions"]["configOverrides"] = overrides
    return decode_canonical_request(
        payload,
        resource_paths={r.resource_id: r.source_path for r in encoded.resources},
        output_directory=tmp_path,
    )


def test_request_decoding_keeps_the_font_leaf_alone(tmp_path):
    # The interval is derived when the override is applied, so a replayed
    # request re-encodes to the same overrides (Gallery replay equivalence).
    decoded = _decode(tmp_path, {FONT: 30})
    overrides = dict(decoded.options.config_overrides)
    assert overrides[FONT] == 30 and INTERVAL not in overrides
    reencoded = encode_canonical_request(decoded).payload["diagramOptions"]
    assert reencoded["configOverrides"] == overrides
    resolved = _resolve_diagram_options_config(decoded.options)
    assert resolved.objects.definition.circular.interval == 32


def test_legacy_flat_definition_font_alias_uses_the_same_rule(tmp_path):
    decoded = _decode(tmp_path, {"circular_definition_font_size": 30})
    overrides = dict(decoded.options.config_overrides)
    assert overrides[FONT] == 30 and INTERVAL not in overrides
    resolved = _resolve_diagram_options_config(decoded.options)
    assert resolved.objects.definition.circular.interval == 32


def _cli_svg(monkeypatch, tmp_path: Path, extra: list[str]) -> str:
    prefix = tmp_path / "cli"
    monkeypatch.setattr(
        sys,
        "argv",
        ["gbdraw", "circular", "--gbk", str(RECORD), "-o", str(prefix), "-f", "svg", *extra],
    )
    try:
        cli.main()
    except SystemExit as exit_:
        assert exit_.code in (0, None)
    return prefix.with_suffix(".svg").read_text(encoding="utf-8")


def _definition(svg: str) -> str:
    start = svg.index('data-gbdraw-role="record-definition"')
    return svg[start : svg.index("</g>", start)]


def test_cli_and_typed_request_render_the_same_definition(monkeypatch, tmp_path):
    cli_svg = _cli_svg(monkeypatch, tmp_path, ["--definition_font_size", "30"])
    encoded = encode_canonical_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(RECORD)),),
            options=CircularDiagramOptions(config_overrides={FONT: 30.0}),
            output=RenderOutputRequest(
                output_prefix="typed", output_directory=tmp_path, formats=("svg",)
            ),
        )
    )
    payload = json.loads(json.dumps(encoded.payload))
    request = decode_canonical_request(
        payload,
        resource_paths={r.resource_id: r.source_path for r in encoded.resources},
        output_directory=tmp_path,
    )
    typed_svg = render_request(request).output_paths[0].read_text(encoding="utf-8")
    assert _definition(typed_svg) == _definition(cli_svg)
    (tmp_path / "default").mkdir()
    default_svg = _cli_svg(monkeypatch, tmp_path / "default", [])
    assert _definition(default_svg) != _definition(cli_svg)
