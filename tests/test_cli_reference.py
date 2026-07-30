from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "docs" / "CLI_Reference.md"
LONG_OPTION = re.compile(r"(?<!\w)--[a-zA-Z0-9_-]+")


def _live_help_options(mode: str) -> set[str]:
    result = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", mode, "--help"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return set(LONG_OPTION.findall(result.stdout))


def _reference_help_options(mode: str) -> set[str]:
    source = REFERENCE.read_text(encoding="utf-8")
    begin = f"<!-- BEGIN GENERATED {mode.upper()} HELP -->"
    end = f"<!-- END GENERATED {mode.upper()} HELP -->"
    assert source.count(begin) == source.count(end) == 1
    block = source.split(begin, 1)[1].split(end, 1)[0]
    return set(LONG_OPTION.findall(block))


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_reference_generated_help_option_sets_match_live_help(mode: str) -> None:
    assert _reference_help_options(mode) == _live_help_options(mode)
