"""``gbdraw auto``: propose a first drawing from the input and draw it."""

from __future__ import annotations

import argparse
import shlex
import sys
from collections.abc import Sequence
from pathlib import Path

from gbdraw.io.genome import load_gbks
from gbdraw.recommend import Proposal, propose


def parse_auto_args(cmd_args: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="gbdraw auto",
        description=(
            "Read GenBank file(s), propose a complete drawing from what they "
            "contain (mode, palette and color rules, feature types, labels, "
            "strands, GC and skew tracks, legend), print the proposal with its "
            "reasons and the equivalent circular or linear command, then draw "
            "it. The tables the proposal needs are written next to the output "
            "as <prefix>.auto.<table>.tsv, so the printed command can be edited "
            "and rerun."
        ),
    )
    parser.add_argument(
        "--gbk",
        nargs="+",
        required=True,
        metavar="GBK_FILE",
        help="GenBank/DDBJ flat file(s)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output file prefix (default: accession number of the first record)",
    )
    parser.add_argument(
        "-f",
        "--format",
        default="svg",
        help="Comma-separated list of output file formats (default: svg)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files (default: refuse to overwrite).",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Draw without asking. Without it, the proposal is drawn only after "
        "a y answer on a terminal; from a script it is printed and not drawn.",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Print the proposal and the equivalent command and write nothing.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the proposal as JSON instead of text.",
    )
    return parser.parse_args(list(cmd_args))


def _mode_command(
    proposal: Proposal,
    gbk_files: Sequence[str],
    *,
    tables: dict[str, Path],
    output_prefix: str,
    formats: str,
    overwrite: bool,
) -> list[str]:
    """The ``gbdraw <mode>`` arguments that draw the proposal."""
    args = ["--gbk", *gbk_files, *proposal.cli_args(tables), "-o", output_prefix, "-f", formats]
    if overwrite:
        args.append("--overwrite")
    return args


def _confirmed(yes: bool) -> bool:
    if yes:
        return True
    if not sys.stdin.isatty():
        print("Not drawn: pass --yes to draw, or --dry_run to only see the proposal.")
        return False
    answer = input("Draw it like this? [Y/n] ").strip().lower()
    return answer in {"", "y", "yes"}


def auto_main(cmd_args: Sequence[str]) -> None:
    """``gbdraw auto``: propose, explain, and draw."""
    args = parse_auto_args(cmd_args)
    records = load_gbks(list(args.gbk))
    proposal = propose(records)
    output_prefix = str(args.output or records[0].id)
    output_dir = Path(output_prefix).parent
    planned_tables = {name: output_dir / f"{Path(output_prefix).name}.auto.{name}.tsv" for name in proposal.tables()}
    mode_args = _mode_command(
        proposal,
        args.gbk,
        tables=planned_tables,
        output_prefix=output_prefix,
        formats=args.format,
        overwrite=args.overwrite,
    )
    command = " ".join(shlex.quote(token) for token in ("gbdraw", proposal.mode, *mode_args))
    if args.json:
        print(proposal.to_json())
    else:
        print(proposal.describe())
        print()
        print(f"Equivalent command:\n  {command}")
        if planned_tables:
            print("  (the .auto.*.tsv tables are written when the drawing is made)")
        print()
    if args.dry_run or not _confirmed(args.yes):
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    proposal.write_tables(output_dir, Path(output_prefix).name)
    if proposal.mode == "circular":
        from gbdraw.circular import circular_main

        circular_main(mode_args)
    else:
        from gbdraw.linear import linear_main

        linear_main(mode_args)


__all__ = ["auto_main", "parse_auto_args"]
