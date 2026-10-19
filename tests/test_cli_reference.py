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


@pytest.mark.parametrize("mode", ["circular", "linear", "render"])
def test_cli_reference_generated_help_option_sets_match_live_help(mode: str) -> None:
    assert _reference_help_options(mode) == _live_help_options(mode)


def test_cli_reference_names_the_default_label_blacklist_keywords_of_their_owners() -> None:
    """OV-470 (Owner Q-02 a): the CLI keeps no default blacklist; the reference
    names the Web's pre-filled keywords and the Python API's config default."""
    import argparse

    from gbdraw.cli_utils.common import add_label_args
    from gbdraw.config.toml import load_config_toml

    reference = " ".join(REFERENCE.read_text(encoding="utf-8").split())
    web = re.search(r"manualBlacklist: '([^']*)'", (ROOT / "gbdraw/web/js/state.js").read_text(encoding="utf-8"))
    assert web and f"starts as `{web.group(1)}`" in reference
    keywords = load_config_toml("gbdraw.data", "config.toml")["labels"]["filtering"]["blacklist_keywords"]
    assert [keyword.strip() for keyword in web.group(1).split(",")] == keywords
    assert f'--label_blacklist "{",".join(keywords)}"' in reference
    assert "`--label_blacklist` has no default in either mode" in reference
    parser = argparse.ArgumentParser()
    add_label_args(parser)
    assert parser.get_default("label_blacklist") == ""
