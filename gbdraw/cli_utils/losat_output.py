"""CLI helpers shared by the Linear and Circular LOSAT options (design 3.2, 3.7)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from gbdraw.core.numeric import parse_ascii_integer
from gbdraw.render.output_paths import commit_staged_outputs, staged_output_files


def parse_positive_int(value: str) -> int:
    """argparse type of the LOSAT thread and genetic-code options."""

    parsed = parse_ascii_integer(str(value).strip())
    if parsed is None or parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def write_losat_output_files(
    output_dir: Path,
    files: Sequence[tuple[str, str]],
    *,
    overwrite: bool,
) -> None:
    """Write ``--losat_output_dir`` files (name, text): all staged, then committed."""

    output_dir.mkdir(parents=True, exist_ok=True)
    targets = [output_dir / name for name, _text in files]
    with staged_output_files(targets) as staged:
        for target, (_name, text) in zip(targets, files):
            with staged[target].open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
        commit_staged_outputs(staged, overwrite=overwrite)


__all__ = ["parse_positive_int", "write_losat_output_files"]
