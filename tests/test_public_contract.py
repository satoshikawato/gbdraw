from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import json
import os
import sys
from pathlib import Path
from unittest.mock import patch

import gbdraw as public_api
import gbdraw.api as library_api
import gbdraw.api.diagram as diagram_api
import gbdraw.api.options as api_options
import gbdraw.api.render as api_render
import gbdraw.circular as circular_cli
import gbdraw.linear as linear_cli


SNAPSHOT_PATH = Path(__file__).parent / "fixtures" / "public_contract.json"


class _ParserCaptured(Exception):
    def __init__(self, parser: argparse.ArgumentParser):
        self.parser = parser


def _json_value(value):
    if value is inspect.Signature.empty:
        return "<empty>"
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    module = getattr(value, "__module__", None)
    name = getattr(value, "__qualname__", getattr(value, "__name__", None))
    if (module, name) == ("pandas.core.frame", "DataFrame"):
        module = "pandas"
    return f"{module}.{name}" if module and name else str(value)


def _hash_json(value) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _capture_parser(module) -> argparse.ArgumentParser:
    def capture(parser, _args=None, _namespace=None):
        raise _ParserCaptured(parser)

    with (
        patch.object(argparse.ArgumentParser, "parse_args", capture),
        patch.object(sys, "argv", ["gbdraw"]),
        patch.dict(os.environ, {"COLUMNS": "80"}),
    ):
        try:
            module._get_args([])
        except _ParserCaptured as captured:
            return captured.parser
    raise AssertionError("CLI parser was not captured")


def _parser_actions(parser: argparse.ArgumentParser) -> list[dict[str, object]]:
    return [
        {
            "action": type(action).__name__,
            "choices": _json_value(action.choices),
            "const": _json_value(action.const),
            "default": _json_value(action.default),
            "dest": action.dest,
            "help": _json_value(action.help),
            "metavar": _json_value(action.metavar),
            "nargs": _json_value(action.nargs),
            "options": list(action.option_strings),
            "required": action.required,
            "type": _json_value(action.type),
        }
        for action in parser._actions
    ]


def _namespace_contract(namespace: argparse.Namespace) -> dict[str, object]:
    return {
        key: {"type": type(value).__name__, "value": _json_value(value)}
        for key, value in sorted(vars(namespace).items())
    }


def _api_contract() -> list[dict[str, object]]:
    contract = []
    for name in public_api.__all__:
        value = getattr(public_api, name)
        try:
            signature = inspect.signature(value)
        except (TypeError, ValueError):
            signature_contract = None
        else:
            signature_contract = {
                "parameters": [
                    {
                        "name": parameter.name,
                        "kind": parameter.kind.name,
                        "default": _json_value(parameter.default),
                        "annotation": _json_value(parameter.annotation),
                    }
                    for parameter in signature.parameters.values()
                ],
                "return": _json_value(signature.return_annotation),
            }
        contract.append(
            {
                "name": name,
                "module": getattr(value, "__module__", None),
                "qualname": getattr(value, "__qualname__", getattr(value, "__name__", None)),
                "signature": signature_contract,
            }
        )
    return contract


def _cli_contract(module, minimal_args: list[str], representative_args: list[str]) -> dict[str, str]:
    parser = _capture_parser(module)
    return {
        "actions_sha256": _hash_json(_parser_actions(parser)),
        "defaults_sha256": _hash_json(_namespace_contract(module._get_args(minimal_args))),
        "representative_sha256": _hash_json(
            _namespace_contract(module._get_args(representative_args))
        ),
    }


def build_contract() -> dict[str, object]:
    return {
        "api": {
            "all": list(public_api.__all__),
            "contract_sha256": _hash_json(_api_contract()),
        },
        "circular_cli": _cli_contract(
            circular_cli,
            ["--gbk", "contract.gb"],
            [
                "--gbk", "contract.gb", "--output", "diagram", "--palette", "orchid",
                "--nt", "AT", "--window", "1000", "--step", "100", "--features",
                "CDS,tRNA", "--labels", "both", "--legend", "bottom", "--format",
                "svg,interactive_svg", "--track_type", "spreadout", "--separate_strands",
            ],
        ),
        "linear_cli": _cli_contract(
            linear_cli,
            ["--gbk", "contract-a.gb"],
            [
                "--gbk", "contract-a.gb", "contract-b.gb", "--blast", "hits.tsv",
                "--output", "diagram", "--palette", "orchid", "--nt", "AT", "--window",
                "1000", "--step", "100", "--features", "CDS,tRNA", "--legend", "bottom",
                "--format", "svg,interactive_svg", "--gc", "--skew",
                "--separate_strands",
            ],
        ),
    }


def test_public_api_and_cli_contract_snapshot() -> None:
    expected = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    assert build_contract() == expected


def test_low_level_api_owners_are_explicit() -> None:
    removed_diagram_reexports = {
        "DEFAULT_SELECTED_FEATURES",
        "assemble_circular_diagram_from_record",
        "assemble_circular_diagram_from_records",
        "assemble_linear_diagram_from_records",
        "build_circular_diagram",
        "build_circular_multi_diagram",
        "build_linear_diagram",
    }
    removed_canvas_and_configurator_reexports = {
        "BlastMatchConfigurator",
        "CircularCanvasConfigurator",
        "DepthConfigurator",
        "FeatureDrawingConfigurator",
        "GcContentConfigurator",
        "GcSkewConfigurator",
        "LegendDrawingConfigurator",
        "LinearCanvasConfigurator",
    }
    removed_convenience_reexports = {
        "DiagramOptions",
        "OutputOptions",
        "TrackOptions",
        "parse_formats",
        "save_figure",
    }
    removed_top_level_names = (
        removed_diagram_reexports
        | removed_canvas_and_configurator_reexports
        | removed_convenience_reexports
    )

    assert removed_top_level_names.isdisjoint(library_api.__all__)
    assert all(not hasattr(library_api, name) for name in removed_top_level_names)
    assert importlib.util.find_spec("gbdraw.api.canvas") is None
    assert importlib.util.find_spec("gbdraw.api.configurators") is None

    assert removed_diagram_reexports <= set(diagram_api.__all__)
    assert api_options.OutputOptions is not None
    assert not hasattr(api_render, "parse_formats")
    assert not hasattr(api_render, "save_figure")

    assert {"draw_circular", "draw_linear", "read_genbank", "read_gff"} <= set(
        public_api.__all__
    )
    assert (
        library_api.CircularDiagramOptions
        is api_options.CircularDiagramOptions
    )
    assert library_api.LinearDiagramOptions is api_options.LinearDiagramOptions
    assert {
        "CircularDiagramOptions",
        "CircularOutputOptions",
        "CircularTrackOptions",
        "LinearDiagramOptions",
        "LinearOutputOptions",
        "LinearTrackOptions",
        "plan_circular_request",
        "plan_linear_request",
    } <= set(library_api.__all__)
    assert library_api.render_to_bytes is api_render.render_to_bytes
    assert library_api.save_figure_to is api_render.save_figure_to
