from __future__ import annotations

import copy
import importlib.util
import inspect
from dataclasses import fields

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw.analysis as analysis
import gbdraw.analysis.collinearity as collinearity
import gbdraw.cli_utils as cli_utils
import gbdraw.cli_utils.common as cli_common
import gbdraw.cli_utils.session as cli_session
import gbdraw.config.modify as config_modify
import gbdraw.diagrams.circular as circular_diagrams
import gbdraw.diagrams.linear as linear_diagrams
import gbdraw.render as render_package
import gbdraw.render.export as render_export
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_circular_diagram_from_records,
    assemble_linear_diagram_from_records,
)
from gbdraw.api.options import DiagramOptions
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.render.interactive_context import build_interactive_svg_context


@pytest.mark.parametrize(
    "name",
    ["add_output_args", "add_stroke_args", "add_legend_args"],
)
def test_unused_cli_helper_exports_are_removed(name: str) -> None:
    assert not hasattr(cli_common, name)
    assert not hasattr(cli_utils, name)
    assert name not in cli_common.__all__
    assert name not in cli_utils.__all__


def test_obsolete_circular_component_shim_is_removed() -> None:
    assert importlib.util.find_spec("gbdraw.circular_diagram_components") is None


def test_obsolete_label_placement_facade_is_removed() -> None:
    assert importlib.util.find_spec("gbdraw.labels.placement") is None


def test_noop_horizontal_label_compatibility_adapter_is_removed() -> None:
    assert importlib.util.find_spec("gbdraw.labels.circular_horizontal") is None


def test_unused_native_collinearity_wrapper_is_removed() -> None:
    assert not hasattr(collinearity, "build_native_collinearity_blocks")
    assert not hasattr(analysis, "build_native_collinearity_blocks")
    assert "build_native_collinearity_blocks" not in collinearity.__all__
    assert "build_native_collinearity_blocks" not in analysis.__all__


def test_unused_session_compatibility_wrapper_is_removed() -> None:
    assert not hasattr(cli_session, "validate_session_override_args")
    assert "validate_session_override_args" not in cli_session.__all__


def test_unused_diagram_save_wrappers_are_removed() -> None:
    assert not hasattr(circular_diagrams, "plot_circular_diagram")
    assert not hasattr(linear_diagrams, "plot_linear_diagram")
    assert "plot_circular_diagram" not in circular_diagrams.__all__
    assert "plot_linear_diagram" not in linear_diagrams.__all__


def test_internal_render_aliases_and_availability_proxy_are_removed() -> None:
    assert not hasattr(render_package, "parse_formats")
    assert not hasattr(render_package, "save_figure")
    assert not hasattr(render_export, "CAIROSVG_AVAILABLE")
    assert not hasattr(cli_common, "CAIROSVG_AVAILABLE")
    assert not hasattr(cli_utils, "CAIROSVG_AVAILABLE")


def test_dead_flat_config_overrides_are_removed() -> None:
    parameters = inspect.signature(modify_config_dict).parameters
    assert {
        "cicular_width_with_labels",
        "blast_color_min",
        "blast_color_max",
    }.isdisjoint(parameters)
    assert not hasattr(config_modify, "suppress_gc_content_and_skew")
    assert "suppress_gc_content_and_skew" not in config_modify.__all__


def test_depth_tick_interval_is_reader_only() -> None:
    current = load_config_toml("gbdraw.data", "config.toml")
    assert "tick_interval" not in current["objects"]["depth"]

    legacy = copy.deepcopy(current)
    legacy["objects"]["depth"].pop("large_tick_interval")
    legacy["objects"]["depth"]["tick_interval"] = 12

    depth = GbdrawConfig.from_dict(legacy).objects.depth
    assert depth.large_tick_interval == pytest.approx(12)
    assert not hasattr(depth, "tick_interval")


def test_feature_table_aliases_are_removed_from_fresh_python_api() -> None:
    option_fields = {item.name for item in fields(DiagramOptions)}
    assert {"feature_table", "feature_table_file"}.isdisjoint(option_fields)

    for function in (
        assemble_circular_diagram_from_record,
        assemble_circular_diagram_from_records,
        assemble_linear_diagram_from_records,
        build_interactive_svg_context,
    ):
        parameters = inspect.signature(function).parameters
        assert {"feature_table", "feature_table_file"}.isdisjoint(parameters)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"label_placement": "on_feature"}, "label_placement"),
        ({"linear_track_layout": "spreadout"}, "linear_track_layout"),
        ({"linear_track_layout": "tuckin"}, "linear_track_layout"),
    ],
)
def test_fresh_config_overrides_reject_removed_linear_values(
    kwargs: dict[str, str],
    message: str,
) -> None:
    current = load_config_toml("gbdraw.data", "config.toml")
    with pytest.raises(ValidationError, match=message):
        modify_config_dict(current, **kwargs)


@pytest.mark.circular
def test_removed_default_keys_remain_tolerated_in_legacy_config() -> None:
    current = load_config_toml("gbdraw.data", "config.toml")
    assert "circular" not in current["objects"]["gc_content"]
    assert "linear" not in current["objects"]["gc_content"]
    assert "circular" not in current["objects"]["gc_skew"]

    legacy = copy.deepcopy(current)
    legacy["objects"]["gc_content"]["circular"] = {"norm_factor": 0.65}
    legacy["objects"]["gc_content"]["linear"] = {"height": 20}
    legacy["objects"]["gc_skew"]["circular"] = {"norm_factor": 0.80}

    assert GbdrawConfig.from_dict(legacy) == GbdrawConfig.from_dict(current)

    record = SeqRecord(
        Seq("ATGC" * 250),
        id="legacy-config",
        name="legacy-config",
    )
    current_svg = assemble_circular_diagram_from_record(
        record,
        config_dict=copy.deepcopy(current),
        legend="none",
    ).tostring()
    legacy_svg = assemble_circular_diagram_from_record(
        record,
        config_dict=legacy,
        legend="none",
    ).tostring()
    assert legacy_svg == current_svg
