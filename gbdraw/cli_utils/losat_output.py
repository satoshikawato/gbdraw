"""CLI helpers shared by the Linear and Circular LOSAT options (design 3.2, 3.7)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from gbdraw.api.options import LOSATN_TASKS
from gbdraw.core.numeric import parse_ascii_integer
from gbdraw.render.output_paths import commit_staged_outputs, staged_output_files


def parse_positive_int(value: str) -> int:
    """argparse type of the LOSAT thread and genetic-code options."""

    parsed = parse_ascii_integer(str(value).strip())
    if parsed is None or parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def add_losatn_task_arg(parser: argparse.ArgumentParser) -> None:
    """``--losatn_task`` of either mode."""

    parser.add_argument(
        '--losatn_task',
        choices=LOSATN_TASKS,
        default=None,
        help=(
            'LOSATN search task with --losat losatn. A native runtime that does not '
            'support the task stops before searching (default: megablast).'
        ),
    )


def add_losat_runtime_args(parser: argparse.ArgumentParser, *, blast_programs: str) -> None:
    """The search runtime of ``--losat`` in either mode: ``--losat_bin`` or ``--ncbi_blast_bin``, and
    ``--losat_threads``; ``blast_programs`` names the BLAST+ programs the mode's ``--losat`` runs."""

    runtime_group = parser.add_mutually_exclusive_group()
    runtime_group.add_argument(
        '--losat_bin',
        metavar='PATH',
        default=None,
        help='Native LOSAT executable for --losat (default: automatic runtime resolution).',
    )
    runtime_group.add_argument(
        '--ncbi_blast_bin',
        metavar='PATH',
        default=None,
        help=(
            f'NCBI BLAST+ executable of the selected --losat program ({blast_programs}) '
            '(default: automatic runtime resolution).'
        ),
    )
    parser.add_argument(
        '--losat_threads',
        metavar='N',
        type=parse_positive_int,
        default=None,
        help='Threads passed to each LOSAT or NCBI BLAST+ job (default: runtime default).',
    )


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
