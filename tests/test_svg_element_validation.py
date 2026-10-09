"""gbdraw builds its SVG elements without svgwrite's attribute validation.

User colors are checked once when they are read (``test_user_colors.py``);
svgwrite's debug mode would check every attribute of every element again.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from svgwrite.validator2 import Full11Validator

from gbdraw import cli

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
_VALIDATOR_METHODS = (
    "check_all_svg_attribute_values",
    "check_svg_attribute_value",
    "check_svg_type",
    "check_valid_children",
)


def _count_validator_calls(monkeypatch: pytest.MonkeyPatch) -> Counter[str]:
    calls: Counter[str] = Counter()
    for name in _VALIDATOR_METHODS:
        original = getattr(Full11Validator, name)

        def counted(self: Full11Validator, *args: Any, _name: str = name, _original: Any = original) -> Any:
            calls[_name] += 1
            return _original(self, *args)

        monkeypatch.setattr(Full11Validator, name, counted)
    return calls


def _run_cli(*args: str) -> None:
    argv = sys.argv
    sys.argv = ["gbdraw", *args]
    try:
        cli.main()
    finally:
        sys.argv = argv


def test_rendering_does_not_run_svgwrite_validation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls = _count_validator_calls(monkeypatch)
    _run_cli(
        "circular",
        "--gbk", str(EXAMPLES / "MellatMJNV.gb"),
        "--labels", "both",
        "-o", str(tmp_path / "circular"),
        "-f", "svg",
    )
    _run_cli(
        "linear",
        "--gbk", str(EXAMPLES / "LvMJNV.gb"), str(EXAMPLES / "TrcuMJNV.gb"),
        "-b", str(EXAMPLES / "LvMJNV.TrcuMJNV.tblastx.out"),
        "--show_labels", "all",
        "--gc", "--skew",
        "-o", str(tmp_path / "linear"),
        "-f", "svg",
    )

    assert (tmp_path / "circular.svg").exists()
    assert (tmp_path / "linear.svg").exists()
    assert calls == Counter()
