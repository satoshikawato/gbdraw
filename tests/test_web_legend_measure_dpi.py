"""Web legend text measurement must use the DPI the renderer resolves."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.api.config import apply_config_overrides
from gbdraw.core.text import calculate_bbox_dimensions


@pytest.fixture(scope="module")
def helpers():
    source = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/app/python-helpers.js").read_text()
    namespace: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


def _measure(helpers, caption, config=None, overrides=None):
    payload = json.loads(
        helpers["measure_legend_text_json"](
            caption, "Arial", 14, json.dumps(config), json.dumps(overrides or {})
        )
    )
    assert "error" not in payload
    return payload["width"]


def test_default_request_measures_at_the_packaged_render_dpi(helpers):
    caption = "Added legend entry"
    render_dpi = apply_config_overrides(None, None).canvas.dpi
    expected, _ = calculate_bbox_dimensions(caption, "Arial", 14, render_dpi)
    assert _measure(helpers, caption) == expected
    assert expected != calculate_bbox_dimensions(caption, "Arial", 14, 72)[0]


@pytest.mark.parametrize("dpi", [72, 150])
def test_request_dpi_override_reaches_the_measurement(helpers, dpi):
    caption = "Added legend entry"
    expected, _ = calculate_bbox_dimensions(caption, "Arial", 14, dpi)
    assert _measure(helpers, caption, overrides={"canvas.dpi": dpi}) == expected


def test_request_config_dpi_reaches_the_measurement(helpers):
    caption = "Added legend entry"
    config = apply_config_overrides(None, {"canvas.dpi": 150})
    from gbdraw.config.modify import config_to_raw_dict

    expected, _ = calculate_bbox_dimensions(caption, "Arial", 14, 150)
    assert _measure(helpers, caption, config=config_to_raw_dict(config)) == expected


def test_helper_without_request_config_uses_the_packaged_render_dpi(helpers):
    caption = "Added legend entry"
    expected, _ = calculate_bbox_dimensions(
        caption, "Arial", 14, apply_config_overrides(None, None).canvas.dpi
    )
    payload = json.loads(helpers["measure_legend_text_json"](caption, "Arial", 14))
    assert payload["width"] == expected
