from __future__ import annotations

import copy
import importlib.util
import inspect
from dataclasses import fields
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.analysis as analysis
import gbdraw.analysis.collinearity as collinearity
import gbdraw.analysis.protein_colinearity as protein_collinearity
import gbdraw.api.diagram as diagram_api
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
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
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


def test_test_only_collinearity_and_legacy_membership_surfaces_are_removed() -> None:
    assert importlib.util.find_spec("gbdraw.io.collinearity") is None

    removed_collinearity_names = {
        "build_orthogroup_collinearity_anchors",
        "collapse_nearby_duplicate_anchors",
        "deduplicate_unit_pair_anchors",
        "protein_hits_to_collinearity_anchors",
        "select_collinearity_anchor_hits",
    }
    assert removed_collinearity_names.isdisjoint(collinearity.__all__)
    assert all(not hasattr(collinearity, name) for name in removed_collinearity_names)

    removed_membership_names = {
        "EvidenceIndex",
        "build_orthogroups_from_protein_hits",
        "build_pair_evidence_index",
        "build_protein_colinearity_comparisons",
        "cap_hits_per_query",
        "expand_orthogroup_membership_from_evidence",
    }
    assert removed_membership_names.isdisjoint(protein_collinearity.__all__)
    assert all(
        not hasattr(protein_collinearity, name)
        for name in removed_membership_names
    )


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


def test_internal_save_figure_warns_before_compatibility_export(
    tmp_path: Path,
) -> None:
    output = tmp_path / "compatibility.svg"
    drawing = Drawing(filename=str(output))

    with pytest.warns(
        DeprecationWarning,
        match=r"save_figure\(\).*gbdraw 0\.16",
    ):
        render_export.save_figure(drawing, ["svg"])

    assert output.exists()


def test_dead_web_style_module_and_color_dialog_aliases_are_removed() -> None:
    web_root = Path(__file__).parents[1] / "gbdraw" / "web"

    assert not (web_root / "js" / "app" / "annotations" / "style-actions.js").exists()
    assert not (web_root / "js" / "app" / "feature-editor" / "color-actions.js").exists()
    for relative_path in (Path("js/state.js"),):
        source = (web_root / relative_path).read_text(encoding="utf-8")
        assert "individualLabel" not in source


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
    option_fields = {
        item.name
        for options_type in (CircularDiagramOptions, LinearDiagramOptions)
        for item in fields(options_type)
    }
    assert {"feature_table", "feature_table_file"}.isdisjoint(option_fields)

    assemblers = (
        assemble_circular_diagram_from_record,
        assemble_circular_diagram_from_records,
        assemble_linear_diagram_from_records,
    )
    for function in assemblers:
        parameters = inspect.signature(function).parameters
        assert {"feature_table", "feature_table_file"}.isdisjoint(parameters)
        assert {"config_dict", "config_overrides"}.isdisjoint(parameters)
        assert parameters["cfg"].annotation == "GbdrawConfig"
        assert parameters["cfg"].default is inspect.Parameter.empty

    context_parameters = inspect.signature(build_interactive_svg_context).parameters
    assert {"feature_table", "feature_table_file"}.isdisjoint(context_parameters)

    assert {
        "build_circular_diagram",
        "build_circular_multi_diagram",
        "build_linear_diagram",
    } <= set(diagram_api.__all__)
    assert {
        "assemble_circular_diagram_from_record",
        "assemble_circular_diagram_from_records",
        "assemble_linear_diagram_from_records",
    }.isdisjoint(diagram_api.__all__)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"labels.linear.placement": "on_feature"}, "labels.linear.placement"),
        (
            {"canvas.linear.track_layout": "spreadout"},
            "canvas.linear.track_layout",
        ),
        (
            {"canvas.linear.track_layout": "tuckin"},
            "canvas.linear.track_layout",
        ),
    ],
)
def test_fresh_config_overrides_reject_removed_linear_values(
    kwargs: dict[str, str],
    message: str,
) -> None:
    current = load_config_toml("gbdraw.data", "config.toml")
    with pytest.raises(ValidationError, match=message):
        modify_config_dict(current, kwargs)


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
        cfg=GbdrawConfig.from_dict(copy.deepcopy(current)),
        legend="none",
    ).tostring()
    legacy_svg = assemble_circular_diagram_from_record(
        record,
        cfg=GbdrawConfig.from_dict(legacy),
        legend="none",
    ).tostring()
    assert legacy_svg == current_svg
