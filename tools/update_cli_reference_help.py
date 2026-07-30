"""Update generated CLI help blocks in docs/CLI_Reference.md."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_PATH = ROOT / "docs" / "CLI_Reference.md"
MODES = ("circular", "linear")


def _marker(mode: str, edge: str) -> str:
    return f"<!-- {edge} GENERATED {mode.upper()} HELP -->"


def _live_help(mode: str) -> str:
    result = subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", mode, "--help"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.replace(
        "usage: cli.py",
        f"usage: gbdraw {mode}",
        1,
    ).rstrip()


def render_reference(source: str) -> str:
    """Return the reference with both generated blocks refreshed."""

    rendered = source
    for mode in MODES:
        begin = _marker(mode, "BEGIN")
        end = _marker(mode, "END")
        if rendered.count(begin) != 1 or rendered.count(end) != 1:
            raise ValueError(f"Expected one generated {mode} help marker pair.")
        prefix, remainder = rendered.split(begin, 1)
        _, suffix = remainder.split(end, 1)
        block = f"{begin}\n\n```text\n{_live_help(mode)}\n```\n\n{end}"
        rendered = f"{prefix}{block}{suffix}"
    return rendered


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Exit unsuccessfully when the generated help blocks are stale.",
    )
    args = parser.parse_args(argv)

    source = REFERENCE_PATH.read_text(encoding="utf-8")
    rendered = render_reference(source)
    if args.check:
        if rendered != source:
            print(
                "docs/CLI_Reference.md has stale generated CLI help; "
                "run tools/update_cli_reference_help.py.",
                file=sys.stderr,
            )
            return 1
        return 0

    REFERENCE_PATH.write_text(rendered, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
