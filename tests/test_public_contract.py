from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import sys
from pathlib import Path
from unittest.mock import patch

import gbdraw.api as public_api
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
        "help_sha256": hashlib.sha256(parser.format_help().encode("utf-8")).hexdigest(),
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
                "svg,interactive-svg", "--track_type", "spreadout", "--separate_strands",
            ],
        ),
        "linear_cli": _cli_contract(
            linear_cli,
            ["--gbk", "contract-a.gb"],
            [
                "--gbk", "contract-a.gb", "contract-b.gb", "--blast", "hits.tsv",
                "--output", "diagram", "--palette", "orchid", "--nt", "AT", "--window",
                "1000", "--step", "100", "--features", "CDS,tRNA", "--legend", "bottom",
                "--format", "svg,interactive-svg", "--show_gc", "--show_skew",
                "--separate_strands",
            ],
        ),
    }


def test_public_api_and_cli_contract_snapshot() -> None:
    expected = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
    assert build_contract() == expected
