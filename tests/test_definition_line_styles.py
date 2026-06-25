from __future__ import annotations

import copy
import re
from types import SimpleNamespace
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.linear as linear_cli_module
import gbdraw.session_io as session_io
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.text import create_text_element
from gbdraw.definition_line_styles import (
    parse_definition_line_style_assignment,
    parse_definition_line_style_overrides,
)
from gbdraw.exceptions import ValidationError
from gbdraw.render.groups.linear import DefinitionGroup


def _record(label: str = "", record_id: str = "record_a") -> SeqRecord:
    record = SeqRecord(Seq("A" * 100), id=record_id)
    if label:
        record.annotations["gbdraw_record_label"] = label
    return record


def _canvas_config() -> SimpleNamespace:
    return SimpleNamespace(length_param="short")


def _linear_definition_config() -> dict:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["show_replicon"] = False
    definition_cfg["show_accession"] = True
    definition_cfg["show_length"] = True
    return config_dict


def _text_with_kind(svg: str, kind: str) -> str:
    match = re.search(
        rf'(<text\b(?=[^>]*data-definition-line-kind="{kind}")[^>]*>.*?</text>)',
        svg,
        re.S,
    )
    assert match is not None
    return match.group(1)


def test_parse_definition_line_style_overrides_valid_and_last_wins() -> None:
    overrides = parse_definition_line_style_overrides(
        [
            "species:weight=normal,color=red,size=12",
            "record_label:font_weight=bold,fill=rgb(1,2,3)",
            "coordinates:weight=700,color=rgba(1,2,3,0.5)",
        ]
    )

    assert overrides == {
        "name": {
            "font_weight": "bold",
            "fill": "rgb(1,2,3)",
            "font_size": 12.0,
        },
        "length": {
            "font_weight": "700",
            "fill": "rgba(1,2,3,0.5)",
        },
    }


def test_parse_definition_line_style_duplicate_properties_are_last_wins() -> None:
    line_key, style = parse_definition_line_style_assignment(
        "name:weight=normal,font_weight=bold,color=red,fill=blue"
    )

    assert line_key == "name"
    assert style == {"font_weight": "bold", "fill": "blue"}


@pytest.mark.parametrize(
    "raw",
    [
        "name",
        "unknown:weight=bold",
        "name:unknown=bold",
        "name:weight=",
        "name:size=0",
        "name:size=nan",
        "name:weight=950",
        "name:color=red,",
    ],
)
def test_parse_definition_line_style_invalid(raw: str) -> None:
    with pytest.raises(ValueError):
        parse_definition_line_style_assignment(raw)


@pytest.mark.parametrize(
    "cmd_args",
    [
        ["--gbk", "dummy.gb", "--definition_line_style", "name"],
        ["--gbk", "dummy.gb", "--definition_line_style", "unknown:weight=bold"],
        ["--gbk", "dummy.gb", "--definition_line_style", "name:weight="],
    ],
)
def test_linear_cli_definition_line_style_validation(cmd_args: list[str]) -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(cmd_args)


def test_linear_cli_definition_line_style_forwards(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    record = _record("Top Label")
    captured: dict[str, Any] = {}
    real_modify_config_dict = linear_cli_module.modify_config_dict

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(
        linear_cli_module,
        "assemble_linear_diagram_from_records",
        lambda *_args, **_kwargs: Drawing(filename=str(tmp_path / "dummy.svg")),
    )

    def fake_modify_config_dict(config_dict, **kwargs):
        captured["line_styles"] = kwargs.get("linear_definition_line_styles")
        return real_modify_config_dict(config_dict, **kwargs)

    monkeypatch.setattr(linear_cli_module, "modify_config_dict", fake_modify_config_dict)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--definition_line_style",
            "name:weight=bold,color=red",
            "--definition_line_style",
            "name:color=rgb(1,2,3)",
            "--definition_line_style",
            "coordinates:size=9",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["line_styles"] == {
        "name": {"font_weight": "bold", "fill": "rgb(1,2,3)"},
        "length": {"font_size": 9.0},
    }


def test_definition_line_style_config_defaults_and_legacy_inheritance() -> None:
    default_cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    default_styles = default_cfg.objects.definition.linear.line_styles

    assert default_styles.name.font_weight == "normal"
    assert default_styles.replicon.font_weight == "normal"
    assert default_styles.accession.font_weight == "normal"
    assert default_styles.length.font_weight == "normal"
    assert default_styles.name.font_size is None
    assert default_styles.name.fill is None

    legacy = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    legacy["objects"]["definition"]["linear"].pop("line_styles", None)
    legacy_cfg = GbdrawConfig.from_dict(legacy)
    assert legacy_cfg.objects.definition.linear.line_styles.name.font_weight is None


def test_definition_line_style_config_validation() -> None:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    config_dict["objects"]["definition"]["linear"]["line_styles"]["name"]["font_size"] = 0

    with pytest.raises(ValidationError):
        GbdrawConfig.from_dict(config_dict)


def test_modify_config_dict_updates_nested_definition_line_styles() -> None:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    config_dict["objects"]["definition"]["linear"]["line_styles"]["name"]["fill"] = "#111111"

    modified = modify_config_dict(
        config_dict,
        linear_definition_line_styles={
            "name": {"font_weight": "bold"},
            "accession": {"font_size": 9.0, "fill": None},
        },
    )
    line_styles = modified["objects"]["definition"]["linear"]["line_styles"]

    assert line_styles["name"]["font_weight"] == "bold"
    assert line_styles["name"]["fill"] == "#111111"
    assert line_styles["accession"]["font_size"] == 9.0
    assert "fill" not in line_styles["accession"]


def test_session_io_definition_line_styles_parse_cli_and_emit_gui_args() -> None:
    parsed = session_io._definition_line_styles_from_cli_args(
        [
            "--definition_line_style",
            "name:weight=bold",
            "--definition-line-style",
            "coordinates:color=rgb(1,2,3)",
        ]
    )
    run_args: list[str] = []
    invocation_args: list[str] = []

    session_io._append_linear_definition_line_style_args(
        run_args,
        invocation_args,
        {
            "linear_definition_line_styles": {
                "name": {"font_weight": "bold", "fill": "#111111"},
                "accession": {"font_weight": "normal"},
                "length": {"font_size": 9},
            }
        },
    )

    assert parsed == {
        "name": {"font_weight": "bold"},
        "length": {"fill": "rgb(1,2,3)"},
    }
    assert run_args == [
        "--definition_line_style",
        "name:weight=bold,color=#111111",
        "--definition_line_style",
        "length:size=9",
    ]
    assert invocation_args == run_args


def test_create_text_element_defaults_and_custom_attrs() -> None:
    default_text = create_text_element("Default", 1, 2, 12, "normal", "Arial").tostring()
    custom_text = create_text_element(
        "Custom",
        1,
        2,
        12,
        "bold",
        "Arial",
        stroke="#111111",
        fill="#222222",
        extra_attrs={"data-definition-line-kind": "name"},
    ).tostring()

    assert 'stroke="none"' in default_text
    assert 'fill="black"' in default_text
    assert 'stroke="#111111"' in custom_text
    assert 'fill="#222222"' in custom_text
    assert 'data-definition-line-kind="name"' in custom_text


def test_definition_group_emits_line_specific_styles_and_data_attrs() -> None:
    config_dict = _linear_definition_config()
    line_styles = config_dict["objects"]["definition"]["linear"]["line_styles"]
    line_styles["name"] = {"font_weight": "bold", "fill": "#111111"}
    line_styles["accession"] = {"fill": "rgb(1,2,3)"}
    line_styles["length"] = {"font_size": 8}

    definition_group = DefinitionGroup(
        _record("Alpha <i>Beta</i>", "AC123"),
        config_dict,
        _canvas_config(),
    )
    svg = definition_group.get_group().tostring()

    name_text = _text_with_kind(svg, "name")
    accession_text = _text_with_kind(svg, "accession")
    length_text = _text_with_kind(svg, "length")

    assert 'font-weight="bold"' in name_text
    assert 'fill="#111111"' in name_text
    assert 'font-style="italic"' in name_text
    assert 'fill="rgb(1,2,3)"' in accession_text
    assert 'font-size="8"' in length_text or 'font-size="8.0"' in length_text


def test_definition_line_specific_font_size_affects_layout(monkeypatch: pytest.MonkeyPatch) -> None:
    config_dict = _linear_definition_config()
    config_dict["objects"]["definition"]["linear"]["show_length"] = False
    line_styles = config_dict["objects"]["definition"]["linear"]["line_styles"]
    line_styles["name"] = {"font_size": 20}
    line_styles["accession"] = {"font_size": 5}
    calls: list[tuple[str, float, str]] = []

    def fake_bbox(text, _font_family, font_size, _dpi, *, font_weight="normal", font_style="normal"):
        calls.append((str(text), float(font_size), str(font_style)))
        return float(font_size) * len(str(text)), float(font_size)

    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_bbox_dimensions",
        fake_bbox,
    )
    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_svg_bbox_dimensions",
        fake_bbox,
    )

    definition_group = DefinitionGroup(_record("AA", "BBBB"), config_dict, _canvas_config())

    assert definition_group.definition_bounding_box_width == pytest.approx(40.0)
    assert definition_group.definition_bounding_box_height == pytest.approx(30.0)
    assert ("AA", 20.0, "normal") in calls
    assert ("BBBB", 5.0, "normal") in calls
