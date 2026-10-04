"""CLI helpers shared by the Linear and Circular LOSAT options (design 3.2, 3.7)."""

from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from gbdraw.render.output_paths import commit_staged_output_file


def parse_positive_int(value: str) -> int:
    """argparse type of the LOSAT thread and genetic-code options."""

    try:
        parsed = int(str(value).strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def write_losat_output_files(
    output_dir: Path,
    files: Sequence[tuple[str, str]],
    *,
    overwrite: bool,
) -> None:
    """Write ``--losat_output_dir`` files (name, text), each staged then committed."""

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, text in files:
        target = output_dir / name
        with TemporaryDirectory(prefix=f".{name}.", dir=target.parent) as temp_name:
            staged_path = Path(temp_name) / name
            with staged_path.open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
            commit_staged_output_file(staged_path, target, overwrite=overwrite)


__all__ = ["parse_positive_int", "write_losat_output_files"]
